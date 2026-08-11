###############################################################################
# extension simulations with shared latent weather, seasonal
# progression, anisotropy, and peak-hour event selection.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))
output_dir <- file.path(project_dir, "output_corrected")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_sim <- as.integer(Sys.getenv("ESH_JOINT_SIM_N", unset = "1000"))
B <- as.integer(Sys.getenv("ESH_JOINT_SIM_B", unset = "499"))
h_factors <- c(0.125, 0.25, 0.5, 1, 2)

reference <- unique(fread(
  file.path(project_dir, "data", "era5_consistent",
            "era5_201506_consistent.csv"),
  select = c("site_id", "lon", "lat")
))[order(site_id)]
n_sites <- nrow(reference)
coords <- project_coordinates_km(reference$lon, reference$lat)
distance <- as.matrix(dist(coords))
graph_operators <- make_graph_operators(reference, h_factors)
graph_names <- sprintf("graph_%d", seq_along(graph_operators))

make_pair_operator <- function(W) {
  diag(W) <- 0
  W <- (W + t(W)) > 0
  storage.mode(W) <- "double"
  list(L = diag(rowSums(W)) - W,
       divisor = 2 * sum(W[upper.tri(W)]))
}

pair_distances <- distance[upper.tri(distance)]
variogram_breaks <- as.numeric(quantile(
  pair_distances, probs = seq(0, 1, length.out = 6), names = FALSE
))
variogram_operators <- lapply(seq_len(5), function(k) {
  lower <- variogram_breaks[k]
  upper <- variogram_breaks[k + 1L]
  in_bin <- if (k == 1L) distance > 0 & distance <= upper else
    distance > lower & distance <= upper
  make_pair_operator(in_bin)
})
variogram_names <- sprintf("variogram_%d", seq_along(variogram_operators))

quadratic <- c(
  setNames(lapply(graph_operators, function(op) {
    list(L = op$L, divisor = 2 * op$weight_sum)
  }), graph_names),
  setNames(variogram_operators, variogram_names)
)
metric_names <- names(quadratic)

covariance_from_distance <- function(d, range_km = 450, nugget = 0.05) {
  (1 - nugget) * exp(-d / range_km) + nugget * diag(n_sites)
}
sigma_isotropic <- covariance_from_distance(distance)
anisotropic_distance <- as.matrix(dist(cbind(coords[, 1] / 800,
                                             coords[, 2] / 300))) * 450
sigma_anisotropic <- covariance_from_distance(anisotropic_distance)
chol_isotropic <- t(chol(sigma_isotropic + 1e-8 * diag(n_sites)))
chol_anisotropic <- t(chol(sigma_anisotropic + 1e-8 * diag(n_sites)))
gradient <- as.numeric(scale(0.8 * coords[, 1] + coords[, 2]))

expected_quadratic <- function(sigma) {
  vapply(quadratic, function(op) {
    sum(diag(op$L %*% sigma)) / op$divisor
  }, numeric(1))
}
gradient_quadratic <- vapply(quadratic, function(op) {
  as.numeric(crossprod(gradient, op$L %*% gradient) / op$divisor)
}, numeric(1))
expected_isotropic <- expected_quadratic(sigma_isotropic)
expected_anisotropic <- expected_quadratic(sigma_anisotropic)

scenarios <- data.table(
  scenario = c("joint_sign_null", "joint_gradient_contraction",
               "seasonal_progression_null", "seasonal_gradient_contraction",
               "anisotropic_null", "anisotropic_amplitude_contraction",
               "peak_selection_null", "peak_selection_contraction"),
  family = rep(c("null", "alternative"), 4),
  design = rep(c("shared latent gradient", "within-month progression",
                 "anisotropic covariance", "peak-hour selection"), each = 2)
)

simulate_ar1 <- function(n, rho) {
  suppressWarnings(as.numeric(arima.sim(
    model = list(ar = rho), n = n, sd = sqrt(1 - rho^2)
  )))
}

classify <- function(mu) {
  q <- quantile(mu, c(0.25, 0.75), names = FALSE, type = 7)
  fifelse(mu >= q[2], "high", fifelse(mu >= q[1], "middle", "low"))
}

simulate_record <- function(n, scenario) {
  seasonal <- grepl("seasonal", scenario)
  peak <- grepl("peak_selection", scenario)
  anisotropic <- grepl("anisotropic", scenario)
  alternative <- grepl("contraction", scenario)

  weather <- simulate_ar1(n, 0.65)
  if (seasonal) {
    mu <- 1.4 * seq(-1, 1, length.out = n) + 0.60 * weather
  } else if (peak) {
    hour <- 0:23
    diurnal <- 1.8 * sin(2 * pi * (hour - 2) / 24)
    hourly_mean <- outer(weather, rep(1, 24)) +
      outer(rep(1, n), diurnal) + matrix(rnorm(n * 24, sd = 0.35), n, 24)
    selected_hour <- max.col(hourly_mean, ties.method = "first")
    mu <- hourly_mean[cbind(seq_len(n), selected_hour)]
  } else {
    mu <- weather
  }
  regime <- classify(mu)

  # The gradient sign shares a latent driver with the regional mean, while its
  # squared energy remains constant under the joint null.
  mu_z <- as.numeric(scale(mu))
  sign_driver <- 0.75 * mu_z + sqrt(1 - 0.75^2) * rnorm(n)
  gradient_sign <- ifelse(sign_driver >= 0, 1, -1)
  gradient_amplitude <- rep(1.8, n)
  anomaly_amplitude <- rep(1, n)
  if (scenario == "joint_gradient_contraction") {
    gradient_amplitude[regime == "high"] <- 1.2
  }
  if (scenario == "seasonal_gradient_contraction") {
    gradient_amplitude <- 1.95 - 0.65 * (seq_len(n) - 1) / (n - 1)
  }
  if (scenario %chin% c("anisotropic_amplitude_contraction",
                        "peak_selection_contraction")) {
    anomaly_amplitude[regime == "high"] <- 0.94
  }

  chol_use <- if (anisotropic) chol_anisotropic else chol_isotropic
  state <- chol_use %*% rnorm(n_sites)
  field <- matrix(NA_real_, n_sites, n)
  for (tt in seq_len(n)) {
    innovation <- chol_use %*% rnorm(n_sites)
    state <- 0.45 * state + sqrt(1 - 0.45^2) * innovation
    anomaly <- state - mean(state)
    field[, tt] <- mu[tt] + anomaly_amplitude[tt] * anomaly +
      gradient_sign[tt] * gradient_amplitude[tt] * gradient
  }

  expected_base <- if (anisotropic) expected_anisotropic else expected_isotropic
  expected <- vapply(seq_len(n), function(tt) {
    anomaly_amplitude[tt]^2 * expected_base +
      gradient_amplitude[tt]^2 * gradient_quadratic
  }, numeric(length(metric_names)))
  expected <- t(expected)
  colnames(expected) <- metric_names
  list(mu = mu, regime = regime, field = field, expected = expected,
       latent_sign_correlation = cor(mu, gradient_sign))
}

field_metrics <- function(field) {
  values <- lapply(quadratic, function(op) {
    colSums(field * (op$L %*% field)) / op$divisor
  })
  setDT(as.data.frame(values, check.names = FALSE))
}

record_effect <- function(values, regime) {
  high <- regime == "high"
  middle <- regime == "middle"
  colMeans(values[high, , drop = FALSE]) /
    colMeans(values[middle, , drop = FALSE]) - 1
}

simulate_dataset <- function(scenario) {
  records <- lapply(seq_len(6), function(r) {
    simulate_record(30 + as.integer(r %% 3L == 2L), scenario)
  })
  observed <- lapply(records, function(x) {
    list(regime = x$regime, outcomes = as.matrix(field_metrics(x$field)),
         expected = x$expected)
  })
  target_by_record <- vapply(records, function(x) {
    record_effect(x$expected, x$regime)
  }, numeric(length(metric_names)))
  target <- rowMeans(target_by_record)
  list(records = observed,
       target_graph = mean(target[graph_names]),
       target_variogram = mean(target[variogram_names]),
       latent_sign_correlation = mean(vapply(
         records, `[[`, numeric(1), "latent_sign_correlation")))
}

joint_shift_test <- function(records, B, seed) {
  observed_records <- vapply(records, function(x) {
    record_effect(x$outcomes, x$regime)
  }, numeric(length(metric_names)))
  observed <- rowMeans(observed_records)
  set.seed(seed)
  shifted <- matrix(0, B, length(metric_names),
                    dimnames = list(NULL, metric_names))
  for (x in records) {
    n <- nrow(x$outcomes)
    high <- x$regime == "high"
    middle <- x$regime == "middle"
    by_offset <- matrix(NA_real_, n, length(metric_names))
    for (offset in 0:(n - 1L)) {
      index <- ((seq_len(n) - 1L + offset) %% n) + 1L
      outcomes <- x$outcomes[index, , drop = FALSE]
      by_offset[offset + 1L, ] <- colMeans(outcomes[high, , drop = FALSE]) /
        colMeans(outcomes[middle, , drop = FALSE]) - 1
    }
    shifted <- shifted + by_offset[sample.int(n, B, replace = TRUE), ]
  }
  shifted <- shifted / length(records)
  graph_observed <- mean(observed[graph_names])
  graph_null <- rowMeans(shifted[, graph_names, drop = FALSE])
  variogram_observed <- mean(observed[variogram_names])
  variogram_null <- rowMeans(shifted[, variogram_names, drop = FALSE])
  list(
    estimate_graph = graph_observed,
    estimate_variogram = variogram_observed,
    p_graph = (1 + sum(graph_null <= graph_observed)) / (B + 1),
    p_variogram = (1 + sum(variogram_null <= variogram_observed)) / (B + 1)
  )
}

set.seed(20260874L)
rows <- vector("list", nrow(scenarios) * n_sim)
counter <- 0L
for (j in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[j]
  cat(sprintf("Joint DGP scenario %s (%d simulations)\n", scenario, n_sim))
  for (s in seq_len(n_sim)) {
    counter <- counter + 1L
    simulated <- simulate_dataset(scenario)
    test <- joint_shift_test(simulated$records, B, 20260874L + counter)
    rows[[counter]] <- cbind(
      scenarios[j],
      data.table(
        simulation = s,
        target_graph = simulated$target_graph,
        estimate_graph = test$estimate_graph,
        graph_error = test$estimate_graph - simulated$target_graph,
        p_graph = test$p_graph,
        target_variogram = simulated$target_variogram,
        estimate_variogram = test$estimate_variogram,
        variogram_error = test$estimate_variogram -
          simulated$target_variogram,
        p_variogram = test$p_variogram,
        latent_sign_correlation = simulated$latent_sign_correlation
      )
    )
  }
}

results <- rbindlist(rows)
summary <- results[, .(
  simulations = .N,
  mean_latent_sign_correlation = mean(latent_sign_correlation),
  target_graph = mean(target_graph),
  bias_graph = mean(graph_error),
  rmse_graph = sqrt(mean(graph_error^2)),
  rejection_graph = mean(p_graph <= 0.05),
  target_variogram = mean(target_variogram),
  bias_variogram = mean(variogram_error),
  rmse_variogram = sqrt(mean(variogram_error^2)),
  rejection_variogram = mean(p_variogram <= 0.05)
), by = .(scenario, family, design)]

fwrite(results, file.path(output_dir, "joint_dgp_simulation_results.csv"))
fwrite(summary, file.path(output_dir, "joint_dgp_simulation_summary.csv"))
print(summary)

###############################################################################
# A short annual-sequence extension reproduces the actual 1991--2025 calendar
# with 2015 and 2022 reserved, rather than compressing the retained summers.
###############################################################################

year_sim_n <- as.integer(Sys.getenv("ESH_IRREGULAR_YEAR_SIM_N", unset = "10000"))
calendar_years <- 1991:2025
retained <- !calendar_years %in% c(2015, 2022)

calendar_lrv <- function(values, years, lag = 2L) {
  centred <- values - mean(values)
  gamma0 <- mean(centred^2)
  lrv <- gamma0
  for (k in seq_len(lag)) {
    left <- match(years + k, years)
    valid <- !is.na(left)
    gamma <- mean(centred[valid] * centred[left[valid]])
    lrv <- lrv + 2 * (1 - k / (lag + 1)) * gamma
  }
  max(lrv, 0)
}

compressed_lrv <- function(values, lag = 2L) {
  centred <- values - mean(values)
  n <- length(values)
  lrv <- mean(centred^2)
  for (k in seq_len(lag)) {
    gamma <- mean(centred[(k + 1):n] * centred[1:(n - k)])
    lrv <- lrv + 2 * (1 - k / (lag + 1)) * gamma
  }
  max(lrv, 0)
}

set.seed(20260875L)
annual_rows <- vector("list", 3 * 2 * year_sim_n)
counter <- 0L
for (rho in c(0, 0.3, 0.6)) {
  for (effect in c(0, -0.07)) {
    for (s in seq_len(year_sim_n)) {
      counter <- counter + 1L
      full <- effect + 0.12 * simulate_ar1(length(calendar_years), rho)
      values <- full[retained]
      years <- calendar_years[retained]
      estimate <- mean(values)
      se_student <- sd(values) / sqrt(length(values))
      se_calendar <- sqrt(calendar_lrv(values, years) / length(values))
      se_compressed <- sqrt(compressed_lrv(values) / length(values))
      critical <- qt(0.975, df = length(values) - 1)
      annual_rows[[counter]] <- data.table(
        rho, effect, simulation = s, estimate,
        cover_student = abs(estimate - effect) <= critical * se_student,
        cover_calendar_lag2 = abs(estimate - effect) <= critical * se_calendar,
        cover_compressed_lag2 = abs(estimate - effect) <=
          critical * se_compressed,
        reject_student = pt(estimate / se_student,
                            df = length(values) - 1) <= 0.05,
        reject_calendar_lag2 = pt(estimate / se_calendar,
                                  df = length(values) - 1) <= 0.05,
        reject_compressed_lag2 = pt(estimate / se_compressed,
                                    df = length(values) - 1) <= 0.05
      )
    }
  }
}
annual <- rbindlist(annual_rows)
annual_summary <- annual[, .(
  simulations = .N,
  bias = mean(estimate - effect),
  rmse = sqrt(mean((estimate - effect)^2)),
  coverage_student = mean(cover_student),
  coverage_calendar_lag2 = mean(cover_calendar_lag2),
  coverage_compressed_lag2 = mean(cover_compressed_lag2),
  rejection_student = mean(reject_student),
  rejection_calendar_lag2 = mean(reject_calendar_lag2),
  rejection_compressed_lag2 = mean(reject_compressed_lag2)
), by = .(rho, effect)]
fwrite(annual_summary,
       file.path(output_dir, "irregular_year_simulation_summary.csv"))
print(annual_summary)
