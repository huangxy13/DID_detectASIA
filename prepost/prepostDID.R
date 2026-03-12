# Preparation ----
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse) 
library(fixest)    
library(magrittr) 
library(parallel)

# setwd('D:/NUS Dropbox/Xiangyuan Huang/github')
setwd('DID_detectASIA/prepost')



# Data simulation ----
dfSim_PrePost <- function(nGroupT = 4, meanC = 3.46, sdC = 0.5, ef = 0.2, 
                          rho = 0.3, sampleSize = 10, 
                          nMonthPre = 12, nMonthPost = 12, nMonthWash = 3,
                          b1 = -0.001, b2 = 0, # b1 is trend, b2 is level shift
                          idSD = 1){
  
  treatment_effect <- -meanC * ef 
  
  group_trends <- expand.grid(
    group = 1:nGroupT,
    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, 0.1), 
           group_lag = purrr::accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           meanC0 = rnorm(1, meanC, sdC)) %>%
    ungroup()
  
  df <- expand.grid(group = 1:nGroupT, 
                    rep = 1:sampleSize,
                    month = 1:(nMonthPre + nMonthWash + nMonthPost)) %>%
    left_join(group_trends, by = c("group", "month")) %>%
    mutate(timeP = ifelse(month > nMonthPre + nMonthWash, 1, 0)) %>% 
    mutate(individual_error = rnorm(n(), mean = 0, sd = idSD),
           mean_val = meanC0 + (b1 * month) + (treatment_effect * timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)
}


# Simulation engine ----
simu_prepost <- function(scenario){
  for (i in 1:nrow(scenario)){
    dfG <- dfSim_PrePost(nGroupT = scenario$nGroupT[i], 
                         meanC = scenario$meanC[i], ef = scenario$ef[i], 
                         rho = scenario$rho[i], sampleSize = scenario$sampleSize[i],
                         idSD = scenario$idSD[i]) %>%
      filter((month <= 12) | (month > 15)) 
    
    fit_full <- tryCatch(
      feols(value ~ month + timeP, data = dfG, cluster = ~group),
      error = function(e) NULL)
    
    if (is.null(fit_full)) return(NULL)
    
    t_orig <- fit_full$coeftable["timeP", "t value"]
    p_val <- fit_full$coeftable["timeP", "Pr(>|t|)"] 
    
    data.frame(scenario[i, ], 
               coef = coef(fit_full)["timeP"], 
               pval = p_val)
  }
}



# Robust pre-post model ----
simu_prepost_robust <- function(scenario, n_sims = 100) {
  final_res <- purrr::map_df(1:nrow(scenario), function(i) {
    sim_results <- replicate(n_sims, {
      dfG <- dfSim_PrePost(
        nGroupT    = scenario$nGroupT[i], 
        meanC      = scenario$meanC[i], 
        ef         = scenario$ef[i], 
        rho        = scenario$rho[i], 
        sampleSize = scenario$sampleSize[i],
        idSD       = scenario$idSD[i]) %>% 
        filter(month <= 12 | month > 15)
      
      # 2. Fit Model
      fit <- tryCatch(
        feols(value ~ month + timeP, data = dfG, cluster = ~group),
        error = function(e) NULL
      )
      
      if (is.null(fit)) return(NA)
      
      as.numeric(fit$coeftable["timeP", "Pr(>|t|)"] <= 0.05)
    })
    
    data.frame(
      scenario[i, ], 
      power = mean(sim_results, na.rm = TRUE) * 100,
      n_sims = n_sims
    )
  })
  
  return(final_res)
}



# Final result ----
simu_prepost_final <- function(scenario, n_sims = 500) {
  purrr::map_df(1:nrow(scenario), function(i) {
    results <- replicate(n_sims, {
      df <- dfSim_PrePost(
        nGroupT    = scenario$nGroupT[i], 
        meanC      = scenario$meanC[i], 
        ef         = scenario$ef[i], 
        rho        = scenario$rho[i], 
        sampleSize = scenario$sampleSize[i],
        idSD       = scenario$idSD[i]
      ) %>% filter(month <= 12 | month > 15)
      
      fit <- tryCatch(
        feols(value ~ month + timeP, data = df, cluster = ~group),
        error = function(e) NULL
      )
      
      if (is.null(fit) || !("timeP" %in% names(coef(fit)))) return(NA)
      
      # Extract p-value
      pval <- fit$coeftable["timeP", "Pr(>|t|)"]
      return(as.numeric(pval <= 0.05))
    })
    
    data.frame(
      scenario[i, ], 
      power = mean(results, na.rm = TRUE) * 100,
      actual_sims = sum(!is.na(results))
    )
  })
}



# Wild bootstrap model ----
simu_prepost_wild <- function(scenario, B = 200) {
 results_grid <- purrr::map_df(1:nrow(scenario), function(i) {
    df <- dfSim_PrePost(
      nGroupT    = scenario$nGroupT[i], 
      meanC      = scenario$meanC[i], 
      ef         = scenario$ef[i], 
      rho        = scenario$rho[i], 
      sampleSize = scenario$sampleSize[i],
      idSD       = scenario$idSD[i]) %>% filter(month <= 12 | month > 15)

    # Basic model
    fit_full <- tryCatch(
      feols(value ~ month + timeP, data = df, cluster = ~group),
      error = function(e) NULL)
    if (is.null(fit_full) || !("timeP" %in% names(coef(fit_full)))) {
      return(data.frame(scenario[i, ], coef = NA, 
                        boot_pval = NA, significant = NA))}
    t_orig <- fit_full$coeftable["timeP", "t value"]
    
    # Null model
    fit_null <- feols(value ~ month, data = df)
    df$resid_null <- residuals(fit_null)
    df$fitted_null <- predict(fit_null)
    
    # 4. Wild Bootstrap (Webb weights for G=4 clusters)
    webb_weights <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
    unique_groups <- unique(df$group)
    
    t_boot <- replicate(B, {
      weights <- sample(webb_weights, length(unique_groups), replace = TRUE)
      names(weights) <- unique_groups
      v_vec <- weights[as.character(df$group)]
      
      df$y_star <- df$fitted_null + df$resid_null * v_vec
      
      f_b <- tryCatch(
        feols(y_star ~ month + timeP, data = df, cluster = ~group),
        error = function(e) NULL)
      if(is.null(f_b)) return(NA)
      return(f_b$coeftable["timeP", "t value"])
    })
    
    # Bootstrap p value
    p_val_wild <- mean(abs(na.omit(t_boot)) >= abs(t_orig))
    
    # Return data
    data.frame(
      scenario[i, ], 
      coef = coef(fit_full)["timeP"], 
      boot_pval = p_val_wild,
      significant = as.numeric(p_val_wild <= 0.05)
    )
  })
  
  return(results_grid)
}



# Scenario ----
scenario_grid <- expand.grid(
  nGroupT = 4,             # Number of clusters
  meanC = 3.46,            # Baseline mean
  sdC = 0.5,
  ef = c(0.2),          # 0% vs 20% effect size
  rho = 0.3,               # Autocorrelation
  sampleSize = c(5, 10, 15, 20), # Participants per cluster
  idSD = c(0.1, 0.5, 1)       # Standard deviation of individual error
)

# Run simulation ----
set.seed(456) 
n_cores <- 100
N_iteration <- 100
t0 <- Sys.time()

final_results <- mclapply(1:N_iteration, function(i) {
  simu_prepost_wild(scenario_grid, B = 300)}, 
  mc.cores = n_cores)

t1 <- Sys.time()
print(t1 - t0)

q()



load('results.Rdata')
# 1. Collapse the list into one dataframe
results_df <- bind_rows(final_results)

# 2. Summarize by Scenario
summary_table <- results_df %>%
  group_by(nGroupT, ef, sampleSize, idSD) %>%
  summarise(
    avg_coef = mean(coef, na.rm = TRUE),
    # This is Power if ef > 0, and Type I Error if ef == 0
    rejection_rate = mean(significant, na.rm = TRUE),
    n_total_sims = n(),
    .groups = 'drop'
  )

print(summary_table)


summary_table %>%
  mutate(Effect = ifelse(ef == 0, "Type I Error (Null)", "Power (0.2 Effect)")) %>%
  ggplot(aes(x = sampleSize, y = rejection_rate, color = as.factor(idSD))) +
  geom_line() + 
  geom_point() +
  facet_wrap(~Effect, scales = "free_y") +
  geom_hline(data = filter(summary_table, ef == 0), 
             aes(yintercept = 0.05), linetype = "dashed", color = "red") +
  labs(title = "Simulation Results: Power and Type I Error",
       x = "Sample Size (per Cluster)",
       y = "Rejection Rate (p < 0.05)",
       color = "Individual SD") +
  theme_minimal()