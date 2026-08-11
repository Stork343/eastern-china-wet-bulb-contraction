###############################################################################
# Finite-sample and asymptotic behaviour of repeated-summer inference.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
output_dir <- file.path(project_dir, "output_corrected")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_sim <- as.integer(Sys.getenv("ESH_YEAR_SIM_N", unset = "10000"))
sample_sizes <- c(20L, 33L, 60L, 120L)
dependence <- c(0, 0.3, 0.6)
innovation_types <- c("gaussian", "skewed", "t3")
effects <- c(null = 0, moderate = -0.035, application = -0.07)
year_sd <- 0.12

row_hac_variance <- function(centered, lag) {
  n <- ncol(centered)
  out <- rowMeans(centered^2)
  if (lag > 0L) {
    for (ell in seq_len(lag)) {
      gamma <- rowMeans(
        centered[, (ell + 1L):n, drop = FALSE] *
          centered[, seq_len(n - ell), drop = FALSE]
      )
      out <- out + 2 * (1 - ell / (lag + 1)) * gamma
    }
  }
  pmax(out, .Machine$double.eps)
}

make_innovations <- function(n_sim, n, type) {
  z <- matrix(rnorm(n_sim * n), n_sim, n)
  if (type == "gaussian") return(z)
  if (type == "skewed") {
    log_sd <- 0.65
    mean_raw <- exp(log_sd^2 / 2)
    sd_raw <- sqrt((exp(log_sd^2) - 1) * exp(log_sd^2))
    return((exp(log_sd * z) - mean_raw) / sd_raw)
  }
  # A common radial multiplier retains cross-year independence within rows
  # only if drawn cellwise; standardisation gives unit innovation variance.
  raw <- matrix(rt(n_sim * n, df = 3), n_sim, n)
  raw / sqrt(3)
}

simulate_cell <- function(n, rho, innovation_type, effect) {
  burn_in <- 200L
  total_n <- n + burn_in
  innovation <- make_innovations(n_sim, total_n, innovation_type)
  values <- matrix(0, n_sim, total_n)
  values[, 1] <- innovation[, 1]
  if (total_n > 1L) {
    for (r in 2:total_n) {
      values[, r] <- rho * values[, r - 1L] +
        sqrt(1 - rho^2) * innovation[, r]
    }
  }
  values <- values[, (burn_in + 1L):total_n, drop = FALSE]
  values <- year_sd * values + effect
  estimates <- rowMeans(values)
  centered <- values - estimates
  student_variance <- rowSums(centered^2) / (n - 1L)
  nw2_variance <- row_hac_variance(centered, min(2L, n - 1L))
  growing_lag <- max(1L, floor(n^(1 / 3)))
  hac_variance <- row_hac_variance(centered, growing_lag)

  student_se <- sqrt(student_variance / n)
  nw2_se <- sqrt(nw2_variance / n)
  hac_se <- sqrt(hac_variance / n)
  critical <- qt(0.975, df = n - 1L)
  p_student <- pt(estimates / student_se, df = n - 1L)
  p_nw2 <- pt(estimates / nw2_se, df = n - 1L)
  p_hac <- pnorm(estimates / hac_se)
  negatives <- rowSums(values < 0)
  p_sign <- pbinom(negatives - 1L, n, 0.5, lower.tail = FALSE)

  true_lrv <- year_sd^2 * (1 + rho) / (1 - rho)
  standardized <- sqrt(n) * (estimates - effect) / sqrt(true_lrv)
  data.table(
    simulations = n_sim,
    bias = mean(estimates - effect),
    rmse = sqrt(mean((estimates - effect)^2)),
    growing_lag = growing_lag,
    coverage_student = mean(abs(estimates - effect) <= critical * student_se),
    coverage_nw2 = mean(abs(estimates - effect) <= critical * nw2_se),
    coverage_hac = mean(abs(estimates - effect) <= qnorm(0.975) * hac_se),
    rejection_student = mean(p_student <= 0.05),
    rejection_nw2 = mean(p_nw2 <= 0.05),
    rejection_hac = mean(p_hac <= 0.05),
    rejection_sign = mean(p_sign <= 0.05),
    rejection_three_025 = mean(
      p_student <= 0.025 & p_nw2 <= 0.025 & p_sign <= 0.025
    ),
    z_mean = mean(standardized),
    z_sd = sd(standardized),
    z_q025 = quantile(standardized, 0.025),
    z_q975 = quantile(standardized, 0.975)
  )
}

set.seed(20260804L)
results <- vector(
  "list",
  length(sample_sizes) * length(dependence) * length(innovation_types) *
    length(effects)
)
counter <- 0L
for (n in sample_sizes) {
  for (rho in dependence) {
    for (innovation_type in innovation_types) {
      for (effect_name in names(effects)) {
        counter <- counter + 1L
        cat(sprintf(
          "R=%d rho=%.1f innovation=%s effect=%s\n",
          n, rho, innovation_type, effect_name
        ))
        results[[counter]] <- cbind(
          data.table(
            sample_size = n,
            rho = rho,
            innovation = innovation_type,
            effect = effect_name,
            true_effect = effects[[effect_name]]
          ),
          simulate_cell(n, rho, innovation_type, effects[[effect_name]])
        )
      }
    }
  }
}

summary <- rbindlist(results)
fwrite(summary, file.path(output_dir, "year_inference_simulation_summary.csv"))
print(summary[sample_size == 33 & effect %chin% c("null", "application")])
