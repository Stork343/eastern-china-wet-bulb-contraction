###############################################################################
# Monte Carlo size check for the corrected classification and circular-shift
# test. The earlier short-series block bootstrap is retained only as a
# documented failed calibration in output_corrected/corrected_null_*.csv.
# Under the null, the mean process used for classification is independent of
# the spatial-variance process used as the outcome.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))

output_dir <- file.path(project_dir, "output_corrected")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

simulate_ar1 <- function(n, rho) {
  as.numeric(arima.sim(model = list(ar = rho), n = n, sd = sqrt(1 - rho^2)))
}

one_null_dataset <- function(lengths, rho_mean = 0.65, rho_outcome = 0.45) {
  lapply(seq_along(lengths), function(i) {
    n <- lengths[i]
    list(
      mean = simulate_ar1(n, rho_mean),
      outcome = simulate_ar1(n, rho_outcome)
    )
  })
}

delta_from_strata <- function(strata) {
  extreme <- moderate <- numeric()
  for (x in strata) {
    q <- quantile(x$mean, c(0.25, 0.75), names = FALSE, type = 7)
    extreme <- c(extreme, x$outcome[x$mean >= q[2]])
    moderate <- c(moderate, x$outcome[x$mean >= q[1] & x$mean < q[2]])
  }
  mean(extreme) - mean(moderate)
}

fast_shift_test <- function(strata, B) {
  observed <- delta_from_strata(strata)
  null <- numeric(B)
  for (b in seq_len(B)) {
    sampled <- lapply(strata, function(x) {
      n <- length(x$mean)
      offset <- sample.int(n, 1L) - 1L
      idx <- ((seq_len(n) - 1L + offset) %% n) + 1L
      list(mean = x$mean[idx], outcome = x$outcome[idx])
    })
    # Restore the original mean series; only the outcome is phase shifted.
    for (i in seq_along(sampled)) sampled[[i]]$mean <- strata[[i]]$mean
    null[b] <- delta_from_strata(sampled)
  }
  list(
    observed = observed,
    p_value = (1 + sum(null <= observed)) / (B + 1)
  )
}

set.seed(20260840L)
n_sim <- 500L
B <- 999L
lengths <- c(30L, 31L, 31L, 30L, 31L, 31L)
results <- vector("list", n_sim)

cat(sprintf("Running %d null simulations with %d random shifts each...\n",
            n_sim, B))
for (s in seq_len(n_sim)) {
  dat <- one_null_dataset(lengths)
  set.seed(20260840L + s)
  fit <- fast_shift_test(dat, B = B)
  results[[s]] <- data.table(
    simulation = s,
    delta = fit$observed,
    p_value = fit$p_value,
    rejects_05 = fit$p_value <= 0.05
  )
  if (s %% 25L == 0L) cat(sprintf("  completed %d/%d\n", s, n_sim))
}

results <- rbindlist(results)
summary <- results[, .(
  simulations = .N,
  rejection_rate_05 = mean(rejects_05),
  mean_delta = mean(delta),
  median_p_value = median(p_value)
)]

fwrite(results, file.path(output_dir, "corrected_shift_calibration.csv"))
fwrite(summary, file.path(output_dir, "corrected_shift_calibration_summary.csv"))
print(summary, digits = 4)
