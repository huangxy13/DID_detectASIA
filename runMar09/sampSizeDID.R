# Preparation ----
rm(list = ls())
library(doParallel)
library(fixest)
library(foreach)
library(magrittr)
library(tidyverse)

set.seed(123)
t0 <- Sys.time()

setwd('DID_detectASIA/runMar09')
# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/DID_detectASIA/runMar09')



# Data generation function ----
dfSim <- function(nGroupT = 2, nGroupC = 2, meanC = 3.46, sdC = 0.5, ef = 0.2, 
                  rho = 0.3, sampleSize = 10, 
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 3,
                  b1 = -0.001, b2 = 0, b3 = 0, b4 = -0.01, b6 = -0.05, b7 = -0.1,
                  idSD = 1, shockSpread = 0.01){
  b5 <- -meanC * ef 
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, shockSpread), group_lag = 
             purrr::accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           # .x is the previous accumulated value, .y is the new group_shock value
           meanC0 = rnorm(1, meanC, sdC)) %>%
    ungroup()
  
  df <- group_trends %>% uncount(sampleSize, .id = 'rep') %>%
    mutate(timeP = ifelse(month > nMonthPre, 1, 0),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)}



# Scenarios setup ----
N_reps <- 1001
scenario <- expand.grid(nGroupT = c(2, 3, 4, 5, 6), 
                        nGroupC = 2, meanC = c(3.46), 
                        sdC = c(0.1),
                        ef = c(0, 0.2), rho = c(0.3),
                        sampleSize = c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30), 
                        idSD = c(0.2, 0.3, 0.4),
                        B = 999,
                        N_reps = 1:N_reps)
cat('Number of scenarios: ', nrow(scenario), '\n')



# Simu engine ----
# Define the simulation logic for ONE iteration of ONE scenario
sim_core <- function(index, scenario) {
  B_val <- scenario$B[index]
  dfRCT <- dfSim(
    nGroupT = scenario$nGroupT[index], 
    nGroupC = scenario$nGroupC[index],
    meanC = scenario$meanC[index],
    sdC = scenario$sdC[index], 
    ef = scenario$ef[index], 
    rho = scenario$rho[index], 
    sampleSize = scenario$sampleSize[index],
    idSD = scenario$idSD[index]) %>% 
    filter((month < 13) | (month > 15)) %>%
    mutate(
      item4 = month * timeP,
      item5 = treat * timeP,        
      item6 = treat * month,      
      item7 = treat * month * timeP)
  
  fit_full <- tryCatch(
    feols(value ~ month + timeP + treat + item4 + item5 + item6 + item7, 
          data = dfRCT, cluster = ~group, nthreads = 1),
    error = function(e) NULL)
  if (is.null(fit_full)) return(NULL)
  t_orig <- fit_full$coeftable["item5", "t value"]
  
  fit_null <- feols(value ~ month + timeP + treat + item4 + item6 + item7, 
                    data = dfRCT, nthreads = 1)
  dfRCT$resid_null <- residuals(fit_null)
  dfRCT$fitted_null <- predict(fit_null)
  
  webb_weights <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
  unique_groups <- unique(dfRCT$group)
  
  t_boot <- replicate(B_val, {
    weights <- sample(webb_weights, length(unique_groups), replace = TRUE)
    names(weights) <- unique_groups
    v_vec <- weights[as.character(dfRCT$group)]
    dfRCT$y_star <- dfRCT$fitted_null + dfRCT$resid_null * v_vec
    
    f_b <- tryCatch(
      feols(y_star ~ month + timeP + treat + item4 + item5 + item6 + item7,
            data = dfRCT, cluster = ~group, nthreads = 1), 
      error = function(e) NULL)
    if(is.null(f_b)) return(NA)
    f_b$coeftable["item5", "t value"]
  })
  
  return(data.frame(
    scenario[index, ], 
    coef = coef(fit_full)["item5"], 
    boot_pval = mean(abs(na.omit(t_boot)) >= abs(t_orig))
  ))
}



# Run with mclapply
t0 <- Sys.time()
results_list <- mclapply(1:nrow(scenario), sim_core,
                         scenario = scenario,
                         mc.cores = detectCores() - 1)
final_results <- bind_rows(results_list, .id = "rep")
save(final_results, file = 'results_final.Rdata')
print(Sys.time() - t0)

q()





# Result loading ----
load('results_final.Rdata')

for(i in 1:ncol(final_results)){
  var <- names(final_results)[i]
  if (var %in% c('rep', 'N_reps', 'coef', 'boot_pval')){
    next} else {
      if(length(unique(final_results[, i])) > 1){
        print(var)
        print(unique(final_results[,i]))
      }
    }
}

# Result summary ----
result <- final_results %>%
  arrange(meanC, nGroupT, nGroupC, sdC, ef, idSD, sampleSize) %>% 
  group_by(meanC, nGroupT, nGroupC, sdC, ef, idSD, sampleSize) %>%
  summarise(nIteration = n(), pos = sum(boot_pval <= 0.05),
            power = 100*pos/nIteration) %>%
  mutate(effectSize = factor(ef),
         idSD = as.character(idSD),
         ef = as.character(ef),
         nGroupT = paste0('No. treat site = ', nGroupT)) %>%
  filter(sdC == 0.1)


ggplot(data = result, aes(x = sampleSize, y = power)) +
  geom_point(aes(col = idSD, shape = ef)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  labs(x = 'Sample Size', y = 'Power (%)', 
       shape = 'Effect Size', col = 'Random error') +
  facet_wrap(~nGroupT)
ggsave('power-SampleSizebynTreat.tiff', dpi = 300, width = 15, height = 12)




ggplot(data = result, aes(x = sampleSize, y = power, group = effectSize)) +
  geom_point(aes(col = effectSize)) +
  # geom_line(aes(col = effectSize)) +
  geom_smooth(aes(col = effectSize)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  labs(col ='Effect Size', x = 'Sample Size', y = 'Power (%)')
ggsave('power-sample.tiff', dpi = 300, width = 10, height = 8)


result <- read.csv('result_manual_wild100.csv') %>%
  arrange(ef, sampleSize) %>% 
  group_by(ef, sampleSize) %>%
  summarise(nIteration = n(), pos = sum(boot_pval <= 0.05),
            power = 100*pos/nIteration) %>%
  mutate(effectSize = factor(ef))



q()
# Delete ----
a <- dfG %>% filter((month < 13)|(month > 15))
ggplot(data = a, aes(x = month, y = value)) +
  geom_point(aes(col = treat)) +
  geom_vline(xintercept = 13) +
  geom_smooth()
geom_smooth(value ~ month)