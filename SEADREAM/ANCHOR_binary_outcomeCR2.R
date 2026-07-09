# Preparation ----
rm(list = ls())
library(brglm2)
library(clubSandwich)
library(doParallel)
library(fixest)
library(foreach)
library(magrittr)
library(parallelly)
library(tidyverse)
library(truncnorm)

set.seed(123)
t0 <- Sys.time()

print(getwd())

cits.model <- "val ~ time + intervention + group + time:intervention + group:intervention + group:time + group:intervention:time"
# intervention: indicator for pre/post intervention period
# group: indicator for intervention/control group
# b1*time refers to baseline trend according to control
# b2*intervention refers to immediate level change in control during intervention period
# b3*group refers baseline level difference between intervention and control sites
# b4*time:intervention refers to change in trend after intervention in control sites
# b5*group:intervention refers to difference in immediate level change between intervention and control sites (DiD for level change)
# b6*group:time refers to difference in baseline trend between intervention and control sites
# b7*group:intervention:month refers to difference in trend change between intervention and control sites (DiD for trend change post-intervention)

### Simulate data
# Function to generate data

# Data generation function ----
dfSim <- function(nGroupT = 6, nGroupC = 3, meanC = 25, sdC = 2, ef = 0.2, ## mean mortality is 25% for BSI, CI width is 2.25, sd is 43
                  rho = 0.3, sampleSize = 10, 
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 2,
                  #b1 = -0.001, b2 = 0, b3 = 0, b4 = -0.01, b6 = -0.05, b7 = -0.1, ## keep at the same from Xiangyuan #
                  b1=-0.001, b2=0, b3=0, b4=-0.01, b6=-0.05, b7=-0.1, ## assign a nonzero value to b4
                  idSD = 3, shockSpread = 0.5){ ## defining shockSpread more heterogenous 
  b5 <- -meanC * ef 
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, shockSpread), group_lag = 
             purrr::accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)), # Account for autocorrelation across time with random noise persisting into next month
           # .x is the previous accumulated value, .y is the new group_shock value
           #meanC0 = rnorm(1, meanC, sdC)) %>%
           meanC0=rtruncnorm(n=1,a = 5, b = 95, meanC, sdC)) %>%
    ungroup()
  
  df <- group_trends %>% uncount(sampleSize, .id = 'rep') %>%
    mutate(timeP = ifelse(month > nMonthPre, 1, 0),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error,
           Avalue = ifelse(value >=0, value, 0),
           #Avalue = 1/(1+exp(-value)),
           event = rbinom(n  = n(), 1,  Avalue/100)) %>%
    ungroup()
  return(df)}

# Scenarios setup ----
N_reps <- 200
scenario <- expand.grid(nGroupT = c(6), 
                        nGroupC = c(3), meanC = c(25), 
                        sdC = c(0.1, 0.5),
                        ef = c(0.2, 0.4), rho = c(0.3),
                        sampleSize = c(10, 14, 18), 
                        idSD = c(0.1,0.2),
                        # B = 999,
                        N_reps = 1:N_reps)
cat('Number of scenarios: ', nrow(scenario), '\n')

# Simu engine ----
# Define the simulation logic for ONE iteration of ONE scenario
sim_core <- function(index, scenario) {
  # B_val <- scenario$B[index]
  dfRCT <- dfSim(
    nGroupT = scenario$nGroupT[index], 
    nGroupC = scenario$nGroupC[index],
    meanC = scenario$meanC[index],
    sdC = scenario$sdC[index], 
    ef = scenario$ef[index], 
    rho = scenario$rho[index], 
    sampleSize = scenario$sampleSize[index],
    idSD = scenario$idSD[index]) %>% 
    filter((month < 13) | (month > 14)) %>% ## month 13, 14 belonging to phase-in
    mutate(item4 = month * timeP,
           item5 = treat * timeP,        
           item6 = treat * month,      
           item7 = treat * month * timeP)
  
  fit <- tryCatch(glm(data = dfRCT, family = binomial('logit'), method = "brglmFit",
             event ~ month + timeP + treat + item4 + item5 + item6 + item7),
             error = function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  res <- tryCatch(coef_test(fit, vcov = 'CR2', cluster = dfRCT$group, 
                   coefs = 'item5', test = 'Satterthwaite'),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  
  data.frame(
    scenario[index, ],
    coef = coef(fit)['item5'],
    p_val = res$p_Satt,
    df = res$df_Satt)
}


cat('Number of cores', detectCores(), '\n')
cat('Number of cores', availableCores(), '\n')


# Run with mclapply
t0 <- Sys.time()
result_final <- data.frame()
for (i in 1:nrow(scenario)){
  cat('Iteration: ', i, '\n')
  result <- sim_core(index = i,
                     scenario = scenario)
  result_final %<>% rbind(result)
  
  save(result_final, file = "DID_detectASIA/SEADREAM/ANCHOR_2_binaryoutcome_results_final.Rdata")
}
print(Sys.time() - t0)

q()