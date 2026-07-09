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

did.model <- "val ~ treat + timeP + treat:timeP"
# timeP:       pre/post intervention period (0 = pre, 1 = post)
# treat:       intervention (1) vs control (0) hospital
# treat:timeP: the difference-in-differences estimand (item5) = DiD level change,
#              the extra pre->post change in the treated arm beyond the control arm.
# NOTE: DiD estimates a single LEVEL shift and assumes parallel trends.
#   Group-specific trends are therefore set to zero in the DGP (b6 = b7 = 0),
#   otherwise the level-only DiD would be biased. Common secular trends (b1, b4)
#   are retained; DiD nets them out because they affect both arms equally.
#   Within-group temporal correlation is kept via the AR(1) group_lag in the DGP
#   and via CR2 clustering on group at the inference stage.

### Simulate data
# Data generation function ----
dfSim <- function(nGroupT = 6, nGroupC = 3, meanC = 25, sdC = 2, ef = 0.2,
                  rho = 0.3, sampleSize = 10,
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 2,
                  b1 = -0.001, b2 = 0, b3 = 0, b4 = -0.01, ## common time effects (cancel in DiD)
                  b6 = 0, b7 = 0,   ## differential trends set to 0 -> parallel trends holds
                  idSD = 3, shockSpread = 0.5){
  b5 <- -meanC * ef   ## the DiD level effect (treat x post)
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    group_by(group) %>% arrange(month) %>%
    ## group_lag: AR(1) shock -> mortality of the same group is correlated across months
    mutate(group_shock = rnorm(n(), 0, shockSpread),
           group_lag = purrr::accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           meanC0 = rtruncnorm(n = 1, a = 5, b = 95, meanC, sdC)) %>%
    ungroup()
  
  df <- group_trends %>% uncount(sampleSize, .id = 'rep') %>%
    mutate(timeP = ifelse(month > nMonthPre, 1, 0),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) +
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) +
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error,
           Avalue = ifelse(value >= 0, value, 0),
           event = rbinom(n = n(), 1, Avalue/100)) %>%
    ungroup()
  return(df)}

# Scenarios setup ----
N_reps <- 200
scenario <- expand.grid(nGroupT = c(6),
                        nGroupC = c(3), meanC = c(25),
                        sdC = c(0.1),
                        ef = c(0.2), rho = c(0.3),
                        sampleSize = c(20),
                        idSD = c(0.2),
                        N_reps = 1:N_reps)
cat('Number of scenarios: ', nrow(scenario), '\n')

# Simu engine ----
sim_core <- function(index, scenario) {
  dfRCT <- dfSim(
    nGroupT = scenario$nGroupT[index],
    nGroupC = scenario$nGroupC[index],
    meanC = scenario$meanC[index],
    sdC = scenario$sdC[index],
    ef = scenario$ef[index],
    rho = scenario$rho[index],
    sampleSize = scenario$sampleSize[index],
    idSD = scenario$idSD[index]) %>%
    filter((month < 13) | (month > 14)) %>%  ## months 13, 14 = phase-in, dropped
    mutate(item5 = treat * timeP)            ## DiD estimand
  
  ## DiD model: level shift only (no trend terms)
  fit <- tryCatch(glm(data = dfRCT, family = binomial('logit'), method = "brglmFit",
                      event ~ treat + timeP + item5),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  ## CR2 clustered on group keeps the within-group serial correlation in the SE
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

# Run ----
t0 <- Sys.time()
result_final <- data.frame()
for (i in 1:nrow(scenario)){
  cat('Iteration: ', i, '\n')
  result <- sim_core(index = i, scenario = scenario)
  result_final %<>% rbind(result)
}
print(Sys.time() - t0)

q()