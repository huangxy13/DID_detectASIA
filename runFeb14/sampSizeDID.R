# Preparation ----
rm(list = ls())
library(doParallel)
library(fixest)
library(foreach)
library(magrittr)
library(tidyverse)

set.seed(123)
t0 <- Sys.time()

setwd('DID_detectASIA/runFeb14')
# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/DID_detectASIA/runFeb14')



# Data generation function ----
dfSim <- function(nGroupT = 2, nGroupC = 2, meanC = 3.46, sdC = 0.5, ef = 0.2, 
                  rho = 0.3, sampleSize = 10, 
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 3,
                  b1 = -0.001, b2 = 0, b3 = -1, b4 = 0, b6 = -0.5, b7 = -0.5,
                  idSD = 1){
  b5 <- -meanC * ef 
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, 0.1), group_lag = 
             purrr::accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           meanC0 = rnorm(1, meanC, sdC)) %>%
    ungroup()
  
  df <- expand.grid(group = 1:(nGroupT + nGroupC), 
                    rep = 1:sampleSize,
                    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    left_join(group_trends, by = c("group", "month")) %>%
    mutate(timeP = ifelse(month > 12, 1, 0),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)}



# Scenarios setup ----
iteration <- 1
B <- 100
scenario <- expand.grid(nGroupT = 2, nGroupC = 2, meanC = 3.46, 
                        sdC = c(0.1, 0.3, 0.5), 
                        iteration = iteration,
                        ef = c(0.2), rho = c(0.3),
                        sampleSize = c(10, 15, 20, 25), 
                        idSD = c(0.05, 0.1))



# Simu engine ----
# Define the simulation logic for ONE iteration of ONE scenario
run_single_sim <- function(i, scenario_row, B) {
  # Generate Data
  dfG <- dfSim(
    nGroupT = scenario_row$nGroupT, 
    nGroupC = scenario_row$nGroupC,
    sdC = scenario_row$sdC, 
    meanC = scenario_row$meanC, 
    ef = scenario_row$ef, 
    rho = scenario_row$rho, 
    sampleSize = scenario_row$sampleSize,
    idSD = scenario_row$idSD
  ) %>% filter((month < 13) | (month > 15))
  
  # Create interaction terms
  dfG_agg <- dfG %>% mutate(
    item1 = month * timeP,
    item2 = timeP * treat,        
    item3 = month * treat,      
    item4 = timeP * treat * month
  )
  
  # Original Fit
  fit_full <- tryCatch(
    feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
          data = dfG_agg, cluster = ~group, nthreads = 1),
    error = function(e) NULL)
  
  if (is.null(fit_full)) return(NULL)
  t_orig <- fit_full$coeftable["item2", "t value"]
  
  # Null Model for Wild Bootstrap
  fit_null <- feols(value ~ month + timeP + treat + item1 + item3 + item4, 
                    data = dfG_agg, nthreads = 1)
  dfG_agg$resid_null <- residuals(fit_null)
  dfG_agg$fitted_null <- predict(fit_null)
  
  # Webb Bootstrap
  webb_weights <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
  unique_groups <- unique(dfG_agg$group)
  
  t_boot <- replicate(B, {
    weights <- sample(webb_weights, length(unique_groups), replace = TRUE)
    names(weights) <- unique_groups
    v_vec <- weights[as.character(dfG_agg$group)]
    dfG_agg$y_star <- dfG_agg$fitted_null + dfG_agg$resid_null * v_vec
    
    f_b <- tryCatch(feols(y_star ~ month + timeP + treat + item1 + item2 + item3 + item4,
                          data = dfG_agg, cluster = ~group, nthreads = 1), 
                    error = function(e) NULL)
    if(is.null(f_b)) return(NA)
    f_b$coeftable["item2", "t value"]
  })
  
  p_val <- mean(abs(na.omit(t_boot)) >= abs(t_orig))
  
  return(data.frame(scenario_row, coef = coef(fit_full)["item2"], boot_pval = p_val))
}

# Run with mclapply correctly
N_reps <- 100 
scenario_expanded <- scenario[rep(seq_len(nrow(scenario)), each = N_reps), ]

# Run Simulation
t0 <- Sys.time()
results_list <- mclapply(1:nrow(scenario_expanded), function(i) {
  # Add a try-catch here to see which iteration fails if any
  res <- try(run_single_sim(i, scenario_expanded[i, ], B = 100))
  if(inherits(res, "try-error")) return(NULL)
  return(res)
}, mc.cores = parallel::detectCores() - 1)

# Combine and Save
final_results <- bind_rows(results_list)
save(final_results, file = 'results_final.Rdata')

print(Sys.time() - t0)

q()
  



# Simulate ----
n_cores <- 100
N_iteration <- 100
results <- mclapply(1:N_iteration, function(i) {
  simu(scenario)}, 
  mc.cores = n_cores)



t1 <- Sys.time()
print(t1 - t0)
q()



# Result load ----
load('results_final.Rdata')

plot_data <- final_results %>%
  group_by(ef, sampleSize, rho, idSD) %>%
  summarise(power = mean(boot_pval <= 0.05) * 100, .groups = "drop") %>%
  mutate(EffectSize = factor(ef))

ggplot(plot_data, aes(x = sampleSize, y = power, color = EffectSize)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
  facet_grid(idSD ~ rho, labeller = label_both) +
  theme_minimal() +
  labs(title = "Power Analysis: Effect of Sample Size and Variance",
       x = "Sample Size", y = "Power (%)")



# Parallel Setup ----
num_cores <- parallel::detectCores()
# options(connections = 2048)
print(num_cores)
# cl <- makeCluster(120)
cl <- makeCluster(10)
registerDoParallel(cl)


final_results <- foreach(
  i = 1:nrow(scenario), .combine = rbind, 
  .packages = c("fixest", "dplyr", 'purrr', 'magrittr')) %dopar% {
  
    dfG <- dfSim(nGroupT = scenario$nGroupT[i], nGroupC = scenario$nGroupC[i],
                 sdC = scenario$sdC[i],
                 meanC = scenario$meanC[i], ef = scenario$ef[i], 
                 rho = scenario$rho[i], sampleSize = scenario$sampleSize[i],
                 idSD = scenario$idSD[i]) %>%
      filter((month < 13)|(month > 15))
    
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
write.csv(final_results, 'result_manual_wildH.csv', row.names = FALSE)

t1 <- Sys.time()
print(t1 - t0)

q()



# Result loading ----
res1 <- read.csv('result_manual_wild500a.csv') %>% mutate(case = 'a')
res2 <- read.csv('result_manual_wild500b.csv') %>% mutate(case = 'b')

result <- rbind(res1, res2)

result <- read.csv('result_manual_wildH.csv')

# Result summary ----
result %<>%
  arrange(ef, idSD, sampleSize) %>% 
  group_by(ef, idSD, sampleSize) %>%
  summarise(nIteration = n(), pos = sum(boot_pval <= 0.05),
            power = 100*pos/nIteration) %>%
  mutate(effectSize = factor(ef))

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