# Preparation ----
rm(list = ls())
library(cowplot)
library(doParallel)
library(fixest)
library(foreach)
library(ggsci)
library(ggplot2)
library(magrittr)
library(tidyverse)

set.seed(123)
t0 <- Sys.time()

setwd('DID_detectASIA/SEADREAM')
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
           month = month - 16,
           mean_val = meanC0 + (b1*month) + (b2*timeP) + (b3*treat) + 
             (b4*month*timeP) + (b5*treat*timeP) + (b6*treat*month) + 
             (b7*treat*month*timeP) + group_lag,
           value = mean_val + individual_error) %>%
    ungroup()
  return(df)}



# Scenarios setup ----
N_reps <- 500
scenario <- expand.grid(nGroupT = c(6), 
                        nGroupC = c(3), meanC = c(120), 
                        sdC = c(36),
                        ef = c(0.2), rho = c(0.3),
                        sampleSize = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), 
                        idSD = c(18),
                        B = 999,
                        N_reps = 1:N_reps)
cat('Number of scenarios: ', nrow(scenario), '\n')



# Simu engine ----
# Define the simulation logic for ONE iteration of ONE scenario
sim_core <- function(index, scenario) {
  # Helper so every return path carries the same columns (needed for bind_rows)
  make_row <- function(coef = NA_real_, boot_pval = NA_real_,
                       n_boot_valid = NA_integer_, status = "ok",
                       error_msg = NA_character_) {
    data.frame(
      scenario[index, ],
      coef = coef,
      boot_pval = boot_pval,
      n_boot_valid = n_boot_valid,
      status = status,
      error_msg = error_msg,
      stringsAsFactors = FALSE)
  }

  # Wrap the whole iteration so ANY unexpected error (e.g. in dfSim) is
  # reported as data instead of being silently swallowed by mclapply.
  tryCatch({
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
      filter((month < -3) | (month >= 0)) %>%
      mutate(
        item4 = month * timeP,
        item5 = treat * timeP,
        item6 = treat * month,
        item7 = treat * month * timeP)

    fit_full <- tryCatch(
      feols(value ~ month + timeP + treat + item4 + item5 + item6 + item7,
            data = dfRCT, cluster = ~group, nthreads = 1),
      error = function(e) e)
    if (inherits(fit_full, "error")) {
      return(make_row(status = "fit_full_failed",
                      error_msg = conditionMessage(fit_full)))
    }
    t_orig <- fit_full$coeftable["item5", "t value"]

    fit_null <- tryCatch(
      feols(value ~ month + timeP + treat + item4 + item6 + item7,
            data = dfRCT, nthreads = 1),
      error = function(e) e)
    if (inherits(fit_null, "error")) {
      return(make_row(coef = coef(fit_full)["item5"],
                      status = "fit_null_failed",
                      error_msg = conditionMessage(fit_null)))
    }
    dfRCT$resid_null <- residuals(fit_null)
    dfRCT$fitted_null <- predict(fit_null)

    webb_weights <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
    unique_groups <- unique(dfRCT$group)

    boot_err <- NULL  # records the last bootstrap error message, if any
    t_boot <- replicate(B_val, {
      weights <- sample(webb_weights, length(unique_groups), replace = TRUE)
      names(weights) <- unique_groups
      v_vec <- weights[as.character(dfRCT$group)]
      dfRCT$y_star <- dfRCT$fitted_null + dfRCT$resid_null * v_vec

      f_b <- tryCatch(
        feols(y_star ~ month + timeP + treat + item4 + item5 + item6 + item7,
              data = dfRCT, cluster = ~group, nthreads = 1),
        error = function(e) e)
      if (inherits(f_b, "error")) {
        boot_err <<- conditionMessage(f_b)
        return(NA)
      }
      f_b$coeftable["item5", "t value"]
    })
    n_boot_valid <- sum(!is.na(t_boot))

    make_row(
      coef = coef(fit_full)["item5"],
      boot_pval = mean(abs(na.omit(t_boot)) >= abs(t_orig)),
      n_boot_valid = n_boot_valid,
      status = if (n_boot_valid < B_val) "boot_partial" else "ok",
      error_msg = if (!is.null(boot_err)) boot_err else NA_character_)
  },
  error = function(e) make_row(status = "unexpected_error",
                               error_msg = conditionMessage(e)))
}



# Run with mclapply
t0 <- Sys.time()
results_list <- mclapply(1:nrow(scenario), sim_core,
                         scenario = scenario,
                         mc.cores = detectCores() - 1)
final_results <- bind_rows(results_list, .id = "rep")

# Report on model failures ----
cat('\n--- Iteration status counts ---\n')
print(table(final_results$status, useNA = 'ifany'))
failed <- final_results %>% filter(status != 'ok')
if (nrow(failed) > 0) {
  cat('\n', nrow(failed), ' of ', nrow(final_results),
      ' iterations had problems. Example messages:\n', sep = '')
  print(head(unique(failed$error_msg[!is.na(failed$error_msg)]), 10))
}

save(final_results, file = 'results_final.Rdata')
print(Sys.time() - t0)

q()





# Result loading ----
load('D:/NUS Dropbox/Xiangyuan Huang/github/DID_detectASIA/SEADREAM/results_final0.RData')
final_results0 <- final_results
load('D:/NUS Dropbox/Xiangyuan Huang/github/DID_detectASIA/SEADREAM/results_final.RData')
final_results %<>% rbind(final_results0)


sum <- final_results %>% group_by(sdC, ef, sampleSize, idSD) %>%
  summarise(nPOS = sum(boot_pval < 0.05),
            sum = n(),
            rate = nPOS/sum) %>% ungroup %>%
  arrange(sdC, idSD, ef, sampleSize) %>%
  mutate(ef = ef*100,
         rate = rate*100,
         ef = factor(ef),
         sdC = factor(sdC),
         idSD = factor(idSD))
labels_sdC  <- c("12" = "Low SD", "24" = "Medium SD", "36" = "High SD")

ggplot(data = sum, aes(x = sampleSize, y = rate, color = ef, group = ef)) +
  geom_line(linewidth = 0.8, aes(col = ef)) +
  geom_point(size = 1.8, aes(shape = ef)) +
  geom_hline(aes(yintercept = 80, linetype = "80% Power"), color = "grey40") +
  geom_hline(aes(yintercept = 5, linetype = "5% Power"), color = "grey40") +
  # geom_hline(yintercept = c(80, 5), linetype = "dashed", color = "grey40") +
  facet_grid(sdC~., labeller = labeller(sdC = labels_sdC)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = c(4, 6, 8, 10, 12, 14, 16), limits = c(4, 16)) +
  scale_color_nejm() +
  labs(
    x = "Sample size (No. patient per month per hospital)",
    y = "Power (%)",
    color = "Effect size (%)",
    shape = 'Effect size (%)',
    linetype = 'Reference lines') +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )+
  theme_minimal_grid() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")
# ggsave('D:/SEADREAMpower.jpg', dpi = 300, height = 10, width = 8)


ggsave(
  'D:/SEADREAMpower.pdf',
  device = "pdf",
  width = 148,
  height = 210,
  units = "mm",
  dpi = 300
)




