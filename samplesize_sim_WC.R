# Sample-size finder for patient-level CITS (absolute effect on raw TTT)
# Dependencies
if(!requireNamespace("lme4")) install.packages("lme4")
if(!requireNamespace("pbapply")) install.packages("pbapply")
if(!requireNamespace("parallel")) install.packages("parallel")
library(lme4); library(pbapply); library(parallel); library(future); library(nlme)

# ------------------------
# 1) Simulate one dataset
# ------------------------
simulate_cits <- function(
    n_months_pre = 12,
    n_months_post = 12,
    n_sites_treated = 6,
    n_sites_control = 2,
    patients_per_month = 30,   # m (same across site-months)
    mean_baseline = 72,        # hours
    sigma_within = 36,         # sd of patient-level residual (hours)
    sigma_month = 6,           # sd of month random intercept (hours)
    did_abs_change = -12,      # absolute change (hours) applied to treated post
    ar1_rho = 0,               # optional AR(1) for month effects; 0 = none
    seed = NULL
) {
  if(!is.null(seed)) set.seed(seed)
  total_months <- n_months_pre + n_months_post
  months <- 1:total_months
  post_flag <- as.integer(months > n_months_pre) # temp code for noting post-intervention
  
  sites <- c(paste0("T", seq_len(n_sites_treated)),
             paste0("C", seq_len(n_sites_control)))
  treated_flag <- c(rep(1, n_sites_treated), rep(0, n_sites_control))
  
  # build rows
  rows <- lapply(seq_along(sites), function(si) {
    site <- sites[si]
    treated <- treated_flag[si]
    # repeat for each month
    do.call(rbind, lapply(months, function(m) {
      data.frame(site = site,
                 site_treated = treated,
                 month = m,
                 post = as.integer(m > n_months_pre),
                 patient = seq_len(patients_per_month),
                 stringsAsFactors = FALSE)
    }))
  })
  dat <- do.call(rbind, rows)
  dat$monthf <- factor(dat$month)
  
  # simulate month random effects (optionally AR1)
  if(ar1_rho == 0) {
    month_re <- rnorm(total_months, mean = 0, sd = sigma_month)
  } else {
    innov_sd <- sigma_month * sqrt(1 - ar1_rho^2)
    month_re <- as.numeric(arima.sim(n = total_months, model = list(ar = ar1_rho), sd = innov_sd))
  }
  dat$month_re <- month_re[dat$month]
  
  # patient-level mean for each row
  # baseline mean, plus month effect; treated post gets did_abs_change
  dat$mu <- mean_baseline + dat$month_re + (dat$site_treated == 1 & dat$post == 1) * did_abs_change
  
  # simulate patient outcomes (normal residual)
  dat$y <- rnorm(nrow(dat), mean = dat$mu, sd = sigma_within)
  
  return(dat)
}


# ------------------------
# 2) Fit and test function
# ------------------------
test_cits <- function(dat, alpha = 0.05) {
  # Mixed model: y ~ site_treated * post + (1 | site:month)
  # We use site:month as a grouping factor for month-level random intercepts that are site-specific
  dat$site_month <- interaction(dat$site, dat$month, drop = TRUE)
  form <- y ~ site_treated * post + month + (1 | site_month) # Coef for site_treated:post is the absolute difference in level change due to intervention
  # m <- try(lmer(form, data = dat, REML = TRUE), silent = TRUE)
  m <- try(lme(y ~ site_treated*post + month, random = ~1|site_month, data = dat, correlation = corAR1(form = ~1|site_month)))
  if(inherits(m, "try-error")) return(NA)
  coefs <- coef(summary(m))
  nm <- grep("site_treated.*:.*post|post.*:.*site_treated", rownames(coefs), value = TRUE)
  if(length(nm) == 0) return(NA)
  # pval <- coefs[nm, "Pr(>|t|)"]
  pval <- coefs[nm, "p-value"]
  return(as.numeric(pval < alpha)) # returns 1 if effect size detected; 0 if not
}

# debugonce(test_cits)
# test_cits(a)

# ------------------------
# 3) Estimate power for given m
# ------------------------
estimate_power_for_m <- function(m, nsim = 500, seed = 1234, ...) {
  
  # Set up parallel computing
  future::plan(multisession, workers = parallel::detectCores() - 1)
  
  res <- furrr::future_map_dbl(seq_len(nsim), function(i) { # sapply with progress bar
    withr::with_seed(i, {
      dat <- simulate_cits(patients_per_month = m, seed = seed + i, ...)
      test_cits(dat)
    })
    }, .progress = TRUE)
  res <- res[!sapply(res, is.na)]
  # res_mat <- do.call(rbind, res)
  mean(res)
}


# ------------------------
# Trial
# ------------------------

a <- simulate_cits(did_abs_change = -1,      
                   ar1_rho = 0.5, 
                   mean_baseline = 3,
                   sigma_within = 1.4,
                   sigma_month = 1,
                   seed = 123)
res <- test_cits(a, alpha = 0.05)
pwr <- estimate_power_for_m(m = 30,
                            # n_sites_control = 2,
                            did_abs_change = -12,      
                            # ar1_rho = 0.5, 
                            mean_baseline = 72,
                            sigma_within = 36,
                            sigma_month = 24,
                            seed = 123)
pwr

# ------------------------
# Create a plot for m per site against power
# ------------------------

df <- data.frame(m = 5:15)
df$pwr <- lapply(df$m, function(x) { estimate_power_for_m(m = x, 
                                                          mean_baseline = 4,
                                                          did_abs_change = -0.3,      
                                                          ar1_rho = 0.5, 
                                                          sigma_within = 2.8,
                                                          sigma_month = 1)})
plot(df$m, df$pwr)
