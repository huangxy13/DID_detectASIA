# Preparation ----
library(boot)
library(doParallel)
library(fixest)
library(foreach)
library(ggplot2)
library(magrittr)
library(tidyverse)

set.seed(123)

# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/')
ifelse(str_detect(getwd(), 'github'),
       setwd('DID_detectASIA'), setwd('DID_detectASIA'))
cat("wd =", getwd(), "\n")



# Data gen function ----
dfSim <- function(nGroupT = 6, nGroupC = 2, eps = 0.1,
                  propC = 0.4, ef = 0.2, slope = 0.001, rho = 0.3,
                  sampleSize = 10, nMonthPre = 12, nMonthPost = 12){
  
  df <- expand.grid(group = 1:(nGroupT + nGroupC)) %>%
    mutate(propC = rnorm(nGroupT + nGroupC, propC, 0.005)) %>%
    # group_by(group) %>%
    tidyr::expand_grid(month = -(nMonthPre - 1):nMonthPost) %>%
    mutate(timeP = ifelse(month <= 0, 'pre', 'post'),
           treat = ifelse(group <= nGroupT, 1, 0)) %>%
    arrange(group, month) %>% group_by(group) %>%
    mutate(logitC0 = qlogis(first(propC)))
  
  df %<>% mutate(
    ar_shock = {
      nT <- n()
      e  <- numeric(nT)
      e[1] <- 0
      if (nT > 1) {
        for (t in 2:nT) {
          e[t] <- rho * e[t - 1] + rnorm(1, 0, eps)
        }
      }
      e
    },
    logit_propC = logitC0 + slope * (month - first(month)) + ar_shock,
    logit_prop = ifelse((treat == 0)|(timeP == 'pre'), 
                        logit_propC, logit_propC - ef),
    prop = plogis(logit_prop)) %>% 
    ungroup()
  
  df %<>% tidyr::expand_grid(rep = 1:sampleSize) %>% 
    mutate(event = rbinom(n(), 1, prob = prop))
  
  return(df)
}

# Restricted Wild Bootstrap Statistic Function ----
# Uses the null_model to ensure the bootstrap distribution represents H0.
boot_did_null_logic <- function(data, indices, null_model) {
  eta_null <- predict(null_model, type = "link")
  u_null   <- data$POS - plogis(eta_null)
  
  clusters <- unique(data$group)
  
  # 2. Mammen Weights (v1, v2), suitable for small number of clusters
  v1 <- (1 - sqrt(5)) / 2
  v2 <- (1 + sqrt(5)) / 2
  p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
  w_g <- sample(c(v1, v2), length(clusters), replace = TRUE, prob = c(p1, 1 - p1))
  names(w_g) <- clusters
  w_i <- w_g[as.character(data$group)]
  
  # 3. Generate synthetic 'Null' outcomes
  eta_star <- eta_null + u_null * w_i
  # We use the binomial draw to maintain the binary nature of the data
  prob_star <- plogis(eta_star)
  pos_star  <- rbinom(n = length(prob_star), size = 1, prob = prob_star)
  
  data_star <- data
  data_star$POS_star <- pos_star
  
  # 4. Refit the FULL model on synthetic Null data
  fit_b <- tryCatch(
    feglm(POS_star ~ D | group + month, data = data_star, family = binomial("logit")),
    error = function(e) NULL
  )
  
  if (is.null(fit_b) || !'D' %in% names(coef(fit_b))) return(NA)
  return(coef(fit_b)["D"])
}

# Simulation parameters setup ----
nGroupT <- 6
nGroupC <- 2
ef <- c(0, 0.3)  # Testing Type I error (Size)
rho <- 0.3
propC <- 0.3
sampleSize <- c(10, 30, 50, 70, 90, 110)
iteration <- 1:150
B <- 300 

scenario <- expand.grid(nGroupT = nGroupT, nGroupC = nGroupC,
                        ef = ef, rho = rho, propC = propC,
                        sampleSize = sampleSize,
                        iteration = iteration)



# Parallel Execution ----
cl <- makeCluster(100) 
registerDoParallel(cl)

final_results <- foreach(
  i = 1:nrow(scenario), .combine = rbind, 
  .packages = c("fixest", "boot", "dplyr", "magrittr")) %dopar% {
    # data preparation
    dfG <- dfSim(nGroupT = scenario$nGroupT[i], nGroupC = scenario$nGroupC[i],
                 ef = scenario$ef[i], rho = scenario$rho[i], propC = scenario$propC[i],
                 sampleSize = scenario$sampleSize[i])
    
    dfG_agg <- dfG %>% mutate(post = ifelse(month > 0, 1, 0), 
                              POS = event == 1, 
                              D = treat * post)
    
    # basic model
    fit_actual <- tryCatch(
      feglm(POS ~ D | group + month, data = dfG_agg, family = binomial("logit")),
      error = function(e) NULL
    )
    if(is.null(fit_actual)) return(data.frame(iteration = i, coef = NA, boot_pval = NA))
    obs_coef <- coef(fit_actual)["D"]
    
    # null model
    fit_null <- feglm(POS ~ 1 | group + month, data = dfG_agg, family = binomial("logit"))
    
    # wild bootstrap
    boot_out <- boot(data = dfG_agg, 
                     statistic = boot_did_null_logic, 
                     R = B, 
                     null_model = fit_null)
    boot_dist <- na.omit(boot_out$t)
    
    # p value
    if(length(boot_dist) > 0) {
      p_val <- sum(abs(boot_dist) >= abs(obs_coef)) / length(boot_dist)
    } else {
      p_val <- NA
    }
    
    return(data.frame(iteration = i, coef = obs_coef, boot_pval = p_val))
  }
stopCluster(cl)



# Final Processing ----
scenario_results <- cbind(scenario, final_results[, c("coef", "boot_pval")])
write.csv(scenario_results, 'did3.csv', row.names = FALSE)

# Evaluation
rejection_rate <- mean(scenario_results$boot_pval < 0.05, na.rm = TRUE)
cat("Empirical Rejection Rate (Size):", rejection_rate * 100, "%\n")
q()



# Result ----
res1 <- read.csv('did3a.csv')
res2 <- read.csv('did3b.csv')
res2 %<>% mutate(iteration = iteration + 100)

res <- res1 %>% rbind(res2)
# res <- read.csv('did3c.csv')

stat <- res %>% group_by(ef, sampleSize) %>%
  summarise(total = n(),
            n = sum(boot_pval < 0.05),
            sig = n/total)

ggplot(data = stat, aes(x = sampleSize, y = sig)) +
  geom_point(aes(group = factor(ef))) +
  geom_line(aes(group = factor(ef), col = factor(ef))) +
  labs(x = 'Sample Size', y = 'Power')
ggsave('sample size.tiff', dpi = 300, width = 10, height = 8)

