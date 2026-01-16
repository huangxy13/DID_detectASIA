# Preparation ----
rm(list = ls())
library(doParallel)
library(fixest)
library(foreach)
library(magrittr)
library(tidyverse)

set.seed(123)
t0 <- Sys.time()

# setwd('D:/NUS Dropbox/Xiangyuan Huang/github')
setwd('DID_detectASIA')



# Data generation function ----
dfSim <- function(nGroupT = 6, nGroupC = 2, meanC = 10, sdC = 1, ef = 0.3, 
                  rho = 0.3, sampleSize = 10, nMonthPre = 12, nMonthPost = 12,
                  b1 = -0.001, b2 = 0.5, b3 = -1, b4 = 0, b6 = -1, b7 = -0.5,
                  idSD = 1){
  b5 <- -meanC * ef 
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = -(nMonthPre - 1):nMonthPost) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, 0.1), group_lag = 
             accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           meanC0 = rnorm(1, meanC, sdC)) %>%
    ungroup()
  
  df <- expand.grid(group = 1:(nGroupT + nGroupC), rep = 1:sampleSize,
                    month = -(nMonthPre - 1):nMonthPost) %>%
    left_join(group_trends, by = c("group", "month")) %>%
    mutate(timeP = ifelse(month > 0, 1, 0),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)}



# Scenarios setup ----
iteration <- 1:500
B <- 500
scenario <- expand.grid(nGroupT = 6, nGroupC = 2, sdC = c(0.2), 
                        iteration = iteration,
                        meanC = c(10), ef = c(0, 0.3), rho = c(0.3),
                        sampleSize = c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15), 
                        idSD = c(3))



# Parallel Setup ----
num_cores <- min(parallel::detectCores() - 2, 10) 
cl <- makeCluster(num_cores)
registerDoParallel(cl)


final_results <- foreach(
  i = 1:nrow(scenario), .combine = rbind, 
  .packages = c("fixest", "dplyr", 'purrr', 'magrittr')) %dopar% {
    
    dfG <- dfSim(nGroupT = scenario$nGroupT[i], nGroupC = scenario$nGroupC[i],
                 sdC = scenario$sdC[i],
                 meanC = scenario$meanC[i], ef = scenario$ef[i], 
                 rho = scenario$rho[i], sampleSize = scenario$sampleSize[i],
                 idSD = scenario$idSD[i])

    dfG_agg <- dfG %>% mutate(
      item1 = month * timeP,
      item2 = timeP * treat,        
      item3 = month * treat,      
      item4 = timeP * treat * month)

    fit_full <- tryCatch(
      feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
            data = dfG_agg, cluster = ~group),
      error = function(e) NULL)
    if (is.null(fit_full) || !("item2" %in% names(coef(fit_full)))) return(NULL)
    t_orig <- fit_full$coeftable["item2", "t value"]
    
    # NULL Model (H0: item2 = 0)
    fit_null <- feols(value ~ month + timeP + treat + item1 + item3 + item4, 
                      data = dfG_agg)
    dfG_agg$resid_null <- residuals(fit_null)
    dfG_agg$fitted_null <- predict(fit_null)
    
    # Webb weights
    webb_weights <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
    unique_groups <- unique(dfG_agg$group)
    
    t_boot <- replicate(
      B, 
      {
        weights <- sample(webb_weights, length(unique_groups), replace = TRUE)
        names(weights) <- unique_groups
        v_vec <- weights[as.character(dfG_agg$group)]
        dfG_agg$y_star <- dfG_agg$fitted_null + dfG_agg$resid_null * v_vec
        
        f_b <- tryCatch(feols(
          data = dfG_agg, cluster = ~group, 
          y_star ~ month + timeP + treat + item1 + item2 + item3 + item4), 
          error = function(e) NULL)
        if(is.null(f_b)) return(NA)
        
        return(f_b$coeftable["item2", "t value"])
      }
      )
  
  # 5. Calculate Bootstrap P-value (Two-tailed)
  p_val <- mean(abs(na.omit(t_boot)) >= abs(t_orig))
  
  data.frame(scenario[i, ], coef = coef(fit_full)["item2"], 
             boot_pval = p_val)
}
stopCluster(cl)



# Summary and Save ----
final_results %<>% arrange(meanC, ef, sampleSize, boot_pval)
write.csv(final_results, 'result_manual_wild.csv', row.names = FALSE)

t1 <- Sys.time()
print(t1 - t0)

q()



# Result summary ----
result <- read.csv('result_manual_wild.csv') %>%
  arrange(ef, sampleSize) %>% 
  group_by(ef, sampleSize) %>%
  summarise(nIteration = n(), pos = sum(boot_pval <= 0.05),
            power = 100*pos/nIteration) %>%
  mutate(effectSize = factor(ef))

ggplot(data = result, aes(x = sampleSize, y = power, group = effectSize)) +
  geom_point(aes(col = effectSize)) +
  geom_line(aes(col = effectSize)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  labs(col ='Effect Size', x = 'Sample Size', y = 'Power (%)')
ggsave('power-sample.tiff', dpi = 300, width = 10, height = 8)


