## Install if needed:
## install.packages("fixest")

library(fixest)

set.seed(123)

#-----------------------------
# 1. Data-generating process
#-----------------------------
# Null effect: true beta_treat = 0
simulate_data_null <- function(G = 30,   # number of clusters
                               T = 10,   # time periods
                               n_i = 20  # units per cluster-time
){
  df <- expand.grid(
    cluster = 1:G,
    time    = 1:T,
    id      = 1:n_i
  )
  
  # Treatment: half of clusters treated after halfway in time
  df$treat <- as.integer(df$cluster <= G/2 & df$time > T/2)
  
  # Cluster and time effects + noise, NO treatment effect
  alpha_g <- rnorm(G, 0, 1)    # cluster RE
  gamma_t <- rnorm(T, 0, 0.2)  # time shocks
  
  df$alpha_g <- alpha_g[df$cluster]
  df$gamma_t <- gamma_t[df$time]
  
  eps <- rnorm(nrow(df), 0, 1)
  
  df$y <- 1 + df$alpha_g + df$gamma_t + eps  # true beta_treat = 0
  
  df
}

#----------------------------------------
# 2. One cluster wild bootstrap p-value
#----------------------------------------
wild_cluster_boot_p <- function(df, B = 999){
  # Full model (with treat)
  fit_full <- feols(
    y ~ treat | factor(cluster) + factor(time),
    data     = df,
    cluster  = ~cluster
  )
  
  beta_hat <- coef(fit_full)[["treat"]]
  se_hat   <- se(fit_full)[["treat"]]
  t_obs    <- beta_hat / se_hat
  
  # Restricted (null) model: no treat term
  fit_null <- feols(
    y ~ 1 | factor(cluster) + factor(time),
    data = df
  )
  
  u_hat  <- resid(fit_null)     # residuals under H0
  y_hat0 <- fitted(fit_null)    # fitted under H0
  
  clusters <- sort(unique(df$cluster))
  G        <- length(clusters)
  
  boot_t <- numeric(B)
  
  for(b in seq_len(B)){
    # Rademacher weights at cluster level
    w_g <- sample(c(-1, 1), G, replace = TRUE)
    names(w_g) <- clusters
    w_i <- w_g[as.character(df$cluster)]
    
    # Bootstrap outcome under null:
    # y* = y_hat_null + residual * w_g
    y_star <- y_hat0 + u_hat * w_i
    df$y_star <- y_star
    
    # Refit full model on bootstrap sample
    fit_b <- feols(
      y_star ~ treat | factor(cluster) + factor(time),
      data    = df,
      cluster = ~cluster
    )
    
    beta_b <- coef(fit_b)[["treat"]]
    se_b   <- se(fit_b)[["treat"]]
    
    boot_t[b] <- beta_b / se_b
  }
  
  # Two-sided wild cluster bootstrap p-value
  p_val <- mean(abs(boot_t) >= abs(t_obs))
  p_val
}

#--------------------------------------------------
# 3. Simulation to check that alpha ≈ 5% under H0
#--------------------------------------------------
set.seed(123)
n_sim  <- 200       # increase to e.g. 1000 for smoother estimate
alpha  <- 0.05
Bboot  <- 999

p_vals <- replicate(n_sim, {
  df_sim <- simulate_data_null()
  wild_cluster_boot_p(df_sim, B = Bboot)
})

# Empirical Type I error at nominal 5%:
mean(p_vals < alpha)
