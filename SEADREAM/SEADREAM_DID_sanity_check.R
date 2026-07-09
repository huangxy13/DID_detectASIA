# Preparation ----
rm(list = ls())
library(brglm2)
library(clubSandwich)
library(magrittr)
library(tidyverse)
library(truncnorm)

set.seed(123)

did.model <- "val ~ treat + timeP + treat:timeP"
# timeP:       pre/post intervention period (0 = pre, 1 = post)
# treat:       intervention (1) vs control (0) hospital
# treat:timeP: the difference-in-differences estimand (item5) = DiD level change,
#              the extra pre->post change in the treated arm beyond the control arm.
# NOTE: DiD estimates a single LEVEL shift and assumes parallel trends
#   (group-specific trend terms b6/b7 are fixed at 0 in the DGP below).
#   This is the same DGP as SEADREAM_DID.R / ANCHOR_binary_outcomeCR2.R, but this
#   script runs ONE simulated trial and prints the fit so the DiD + CR2/Satterthwaite
#   pipeline can be sanity-checked before scaling up to a full scenario sweep.

# Data generation function ----
dfSim <- function(nGroupT = 6, nGroupC = 3, meanC = 25, sdC = 2, ef = 0.2,
                  rho = 0.3, sampleSize = 10,
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 2,
                  b1 = -0.001, b2 = 0, b3 = 0, b4 = -0.01, ## common time effects (cancel in DiD)
                  b6 = 0, b7 = 0,   ## differential trends set to 0 -> parallel trends holds
                  idSD = 3, shockSpread = 0.5){
  b5 <- -ef   ## the DiD level effect (treat x post)
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

# Single-scenario parameters ----
## same defaults used across the scenario grid in SEADREAM_DID.R
nGroupT    <- 5
nGroupC    <- 4
meanC      <- 25
sdC        <- 1
ef         <- 10
rho        <- 0.3
sampleSize <- 15
idSD       <- 1

# Simulate one trial ----
dfRCT <- dfSim(nGroupT = nGroupT, nGroupC = nGroupC, meanC = meanC, sdC = sdC,
               ef = ef, rho = rho, sampleSize = sampleSize, idSD = idSD) %>%
  filter((month < 13) | (month > 14)) %>%  ## months 13, 14 = phase-in, dropped
  mutate(item5 = treat * timeP)            ## DiD estimand

cat('Number of clusters (groups):', length(unique(dfRCT$group)), '\n')
cat('Number of observations:', nrow(dfRCT), '\n\n')

# Fit DiD model ----
fit <- glm(data = dfRCT, family = binomial('logit'), method = "brglmFit",
           event ~ treat + timeP + item5)

cat('--- Model coefficients (log-odds scale) ---\n')
print(summary(fit)$coefficients)

# CR2 clustered on group + Satterthwaite df (small-cluster adjustment) ----
res <- coef_test(fit, vcov = 'CR2', cluster = dfRCT$group,
                  coefs = 'item5', test = 'Satterthwaite')

cat('\n--- DiD estimate (item5 = treat:timeP) with CR2/Satterthwaite ---\n')
print(res)

cat('\nInterpretation: coefficient is the log-odds DiD effect; res$df_Satt is the\n')
cat('small-sample-adjusted degrees of freedom (<< number of clusters - 1 typically);\n')
cat('res$p_Satt is the corresponding p-value for the DiD term.\n')
