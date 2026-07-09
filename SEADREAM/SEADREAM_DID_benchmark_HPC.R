# Preparation ----
rm(list = ls())
library(brglm2)
library(clubSandwich)
library(magrittr)
library(parallel)
library(parallelly)
library(tidyverse)
library(truncnorm)

## Force single-threaded BLAS/OpenMP BEFORE forking -- on HPC nodes R is often
## linked against a multi-threaded BLAS (OpenBLAS/MKL). Without this, each of
## the n_cores forked workers can itself spawn several math threads for the
## CR2 matrix computations, so n_cores processes x several threads each fight
## over the same physical cores (oversubscription), which can make parallel
## runs as slow as, or slower than, serial. Setting this here (in the parent,
## before mclapply forks) is inherited by every forked child.
if (requireNamespace('RhpcBLASctl', quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

## L'Ecuyer-CMRG gives each forked worker its own valid, reproducible RNG
## stream (required for correct parallel RNG with mclapply -- see ?mclapply)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
t0 <- Sys.time()

# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/')
setwd('DID_detectASIA/SEADREAM')
print(getwd())

did.model <- "val ~ treat + timeP + treat:timeP"
# timeP:       pre/post intervention period (0 = pre, 1 = post)
# treat:       intervention (1) vs control (0) hospital
# treat:timeP: the difference-in-differences estimand (item5) = DiD level change,
#              the extra pre->post change in the treated arm beyond the control arm.
# NOTE: DiD estimates a single LEVEL shift and assumes parallel trends
#   (group-specific trend terms b6/b7 are fixed at 0 in the DGP below).
#
# This script promotes the single-trial check in SEADREAM_DID_sanity_check.R to a
# 500-iteration simulation, parallelized across HPC cores, using the parameter
# set validated there as the BENCHMARK scenario (nGroupT = 5, nGroupC = 4,
# ef = 10, sampleSize = 15, idSD = 1, sdC = 1). Keep these values fixed here --
# this run is meant as the reference point other scenarios get compared against.

# Data generation function ----
dfSim <- function(nGroupT = 6, nGroupC = 3, meanC = 25, sdC = 2, ef = 0.2,
                  rho = 0.3, sampleSize = 10,
                  nMonthPre = 12, nMonthPost = 12, nMonthWash = 2,
                  b1 = -0.001, b2 = 0, b3 = 0, b4 = -0.01, ## common time effects (cancel in DiD)
                  b6 = 0, b7 = 0,   ## differential trends set to 0 -> parallel trends holds
                  idSD = 3, shockSpread = 0.1){
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

# Benchmark scenario setup ----
## fixed parameters validated in SEADREAM_DID_sanity_check.R -- this is the
## reference scenario, not a sweep, so every row below is identical apart
## from the replicate index
N_reps <- 200
scenario <- expand.grid(nGroupT = 6,
                        nGroupC = 3, meanC = 25,
                        sdC = c(1.5, 2),
                        ef = c(5, 8, 10), rho = 0.3,
                        sampleSize = 25,
                        idSD = 0.5,
                        N_reps = 1:N_reps)
cat('Number of scenarios (iterations): ', nrow(scenario), '\n')

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
  
  data.frame(scenario[index, ],
             coef = coef(fit)['item5'],
             p_val = res$p_Satt,
             df = res$df_Satt)}

# Parallel setup (HPC) ----
## parallelly::availableCores() reads the scheduler's allocation (e.g. Slurm
## SLURM_CPUS_PER_TASK / cgroup limits) instead of the full physical core count
n_cores <- availableCores()
cat('Number of cores available: ', n_cores, '\n')

## mclapply forks the current, already-loaded R process (fast, no socket
## communication, no re-loading packages per worker) instead of spawning new
## Rscript sessions like PSOCK/doSNOW do -- this is the fix for HPC nodes where
## per-worker package reloads over a networked/shared filesystem make PSOCK
## painfully slow to even start.
## NOTE: mclapply is Unix/Linux only -- passing mc.cores > 1 on Windows raises
## an error (rather than silently running serially), so force mc.cores = 1
## there for local syntax-checking; only the HPC (Linux) run is actually parallel.
mc_cores <- if (.Platform$OS.type == "windows") 1 else n_cores

## split the scenario rows into one chunk per core; each chunk is handed to a
## single forked worker, which walks its share with a plain for loop instead
## of relying on mclapply's internal mc.preschedule batching
chunks <- parallel::splitIndices(nrow(scenario), n_cores)
cat('Number of chunks: ', length(chunks), '\n')

# Run ----
t0 <- Sys.time()
result_list <- mclapply(chunks, function(idx) {
  chunk_result <- data.frame()
  for (i in idx) {
    cat('Executing: ', i, '\n')
    result <- sim_core(i, scenario)
    chunk_result %<>% rbind(result)
  }
  chunk_result
}, mc.cores = mc_cores, mc.set.seed = TRUE)
result_final <- do.call(rbind, result_list)
print(Sys.time() - t0)

# Save ----
save(result_final, file = "SEADREAM_DID_benchmark_results.Rdata")

sum <- result_final %>% group_by(meanC, sdC, idSD, ef, rho, sampleSize) %>% 
  summarise(nPOS = sum(p_val < 0.05), nTotal = n(), rate = 100*nPOS/nTotal)
print(sum)


cat('\n--- Benchmark summary ---\n')
cat('Iterations completed: ', nrow(result_final), ' / ', N_reps, '\n')
cat('Power (share with p_val <= 0.05): ',
    mean(result_final$p_val <= 0.05, na.rm = TRUE) * 100, '%\n')
cat('Mean Satterthwaite df: ', mean(result_final$df, na.rm = TRUE), '\n')

q()
