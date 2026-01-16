# Preparation ----
rm(list = ls())
library(doParallel)
library(fixest)
library(foreach)
library(magrittr)
library(tidyverse)

set.seed(123)
t0 <- Sys.time()
print('latest')

# Data gen function (Unchanged logic) ----
dfSim <- function(nGroupT = 6, nGroupC = 2, meanC = 10, sdC = 3, ef = 0.3, 
                  rho = 0.3, sampleSize = 10, nMonthPre = 12, nMonthPost = 12,
                  b1 = -0.001, b2 = 0.5, b3 = -1, b4 = 0, b6 = -1, b7 = -0.5){
  b5 <- -meanC * ef 
  group_trends <- expand.grid(
    group = 1:(nGroupT + nGroupC),
    month = -(nMonthPre - 1):nMonthPost) %>%
    group_by(group) %>% arrange(month) %>%
    mutate(group_shock = rnorm(n(), 0, 1), group_lag = 
             accumulate(group_shock, ~ .x * rho + .y * sqrt(1 - rho^2)),
           meanC0 = rnorm(1, meanC, 0.5)) %>%
    ungroup()
  
  df <- expand.grid(group = 1:(nGroupT + nGroupC),
                    month = -(nMonthPre - 1):nMonthPost,
                    rep = 1:sampleSize) %>%
    left_join(group_trends, by = c("group", "month")) %>%
    mutate(timeP = ifelse(month > 0, 1, 0), # Simplified check
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    mutate(individual_error = rnorm(n(), mean = 0, sd = 0.5),
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)
}

# Restricted Wild Bootstrap Function ----
boot_did_linear <- function(data, null_model) {
  resid_null <- residuals(null_model)
  fit_null <- predict(null_model)
  clusters <- unique(data$group)
  
  v1 <- (1 - sqrt(5)) / 2
  v2 <- (1 + sqrt(5)) / 2
  p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
  
  w_g <- sample(c(v1, v2), length(clusters), replace = TRUE, prob = c(p1, 1 - p1))
  names(w_g) <- clusters
  w_i <- w_g[as.character(data$group)]
  
  data$value_star <- fit_null + resid_null * w_i
  
  # Ensure all items are present in data passed here
  fit_b <- tryCatch(
    feols(value_star ~ month + timeP + treat + item1 + item2 + item3 + item4,  
          data = data, cluster = ~group),
    error = function(e) NULL
  )
  
  if (is.null(fit_b) || !("item2" %in% names(coef(fit_b)))) return(NA)
  
  return(coef(fit_b)["item2"] / se(fit_b)["item2"])
}


boot_naive_cluster <- function(data) {
  clusters <- unique(data$group)
  
  # 1. Resample cluster IDs with replacement
  resampled_ids <- sample(clusters, length(clusters), replace = TRUE)
  
  # 2. Reconstruct data (must create unique group IDs for duplicates)
  data_boot <- lapply(seq_along(resampled_ids), function(i) {
    cluster_df <- data[data$group == resampled_ids[i], ]
    cluster_df$new_group_id <- i  # Treat each selection as a unique cluster
    return(cluster_df)
  }) %>% bind_rows()
  
  # 3. Fit the model on resampled data
  fit_b <- tryCatch(
    feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
          data = data_boot),
    error = function(e) NULL
  )
  
  if (is.null(fit_b) || !("item2" %in% names(coef(fit_b)))) return(NA)
  
  return(coef(fit_b)["item2"])
}



# Parameters ----
iteration <- 1:100
B <- 250 
scenario <- expand.grid(nGroupT = 6, nGroupC = 2, meanC = 10,
                        ef = c(0, 0.3), rho = 0.3,
                        sampleSize = c(10, 15, 20),
                        iteration = iteration)

# Parallel Setup ----
num_cores <- min(parallel::detectCores() - 2, 100)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("dfSim", "boot_naive_cluster"))

# Stage 3: Execution ----
final_results <- foreach(
  i = 1:nrow(scenario), .combine = rbind, 
  .packages = c("fixest", "dplyr", "purrr", "magrittr")) %dopar% {
    
    # Generate data
    dfG <- dfSim(nGroupT = scenario$nGroupT[i], nGroupC = scenario$nGroupC[i],
                 meanC = scenario$meanC[i], ef = scenario$ef[i], 
                 rho = scenario$rho[i], sampleSize = scenario$sampleSize[i])
    
    # Construct model terms
    dfG_agg <- dfG %>% mutate(
      item1 = month * timeP,
      item2 = timeP * treat,        
      item3 = month * treat,      
      item4 = timeP * treat * month)
    
    # # 1. ACTUAL MODEL
    # fit_actual <- tryCatch(
    #   feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
    #         data = dfG_agg, cluster = ~group),
    #   error = function(e) NULL)
    # 
    # if (is.null(fit_actual) || !("item2" %in% names(coef(fit_actual)))) {
    #   return(data.frame(scenario[i, ], coef = NA, boot_pval = NA))
    # }
    # 
    # t_orig <- coef(fit_actual)["item2"] / se(fit_actual)["item2"]
    # 
    # # 2. NULL MODEL (Excluding item2)
    # fit_null <- tryCatch(
    #   feols(value ~ month + timeP + treat + item1 + item3 + item4, data = dfG_agg),
    #   error = function(e) NULL)
    # 
    # if (is.null(fit_null)) return(data.frame(scenario[i, ], coef = NA, boot_pval = NA))
    # 
    # # 3. BOOTSTRAP
    # boot_t_stats <- replicate(B, boot_did_linear(dfG_agg, fit_null))
    # boot_t_stats <- na.omit(boot_t_stats)
    # 
    # p_val <- if(length(boot_t_stats) > 0) {
    #   mean(abs(boot_t_stats) >= abs(t_orig)) 
    # } else { NA }
    # 
    # data.frame(scenario[i, ],
    #            coef = coef(fit_actual)["item2"], 
    #            boot_pval = p_val)
    
    
    # Fit original model to get the observed coefficient
    fit_actual <- tryCatch(
      feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
            data = dfG_agg, cluster = ~group),
      error = function(e) NULL)
    
    if (is.null(fit_actual)) return(NULL)
    obs_coef <- coef(fit_actual)["item2"]
    
    # Run Naive Bootstrap (No Null Model needed)
    boot_coefs <- replicate(B, boot_naive_cluster(dfG_agg))
    boot_coefs <- na.omit(boot_coefs)
    
    # Calculate p-value based on the distribution of coefficients
    # This checks how much of the distribution is on the opposite side of 0
    p_val <- if(length(boot_coefs) > 0) {
      2 * min(mean(boot_coefs <= 0), mean(boot_coefs >= 0))
    } else { NA }
    
    data.frame(scenario[i, ],
               coef = obs_coef, 
               boot_pval = p_val)
  }

stopCluster(cl)

# Cleanup and Save ----
write.csv(final_results, 'DID_detectASIA/result.csv', row.names = FALSE)
t1 <- Sys.time()
print(t1 - t0)



# Result summary ---
final_results %>% group_by(ef, sampleSize) %>%
  summarise(n = n(),
            count = sum(boot_pval <= 0.05))
print(final_results)



q()

# Evaluation
rejection_rate <- mean(scenario_results$boot_pval < 0.05, na.rm = TRUE)
cat("Empirical Rejection Rate (Size):", rejection_rate * 100, "%\n")
q()



# Result ----




# Delete ----
dfG_agg <- df %>% mutate(
  item1 = month * timeP,
  item2 = timeP * treat,  
  item3 = month * treat,      
  item4 = timeP * treat * month)

m <- feols(value ~ month + timeP + treat + item1 + item2 + item3 + item4, 
           data = dfG_agg, cluster = ~group)
summary(m)

m <- feols(value ~ month + timeP + treat, 
           data = dfG_agg, cluster = ~group)
summary(m)