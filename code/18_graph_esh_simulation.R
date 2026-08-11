###############################################################################
# Comprehensive finite-sample simulation for multi-scale graph dispersion.
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

h_factors <- c(0.125, 0.25, 0.5, 1, 2)
graph_names <- sprintf("graph_h_%s", format(h_factors, trim = TRUE))
n_sim <- as.integer(Sys.getenv("ESH_FIELD_SIM_N", unset = "1000"))
B <- as.integer(Sys.getenv("ESH_FIELD_SIM_B", unset = "499"))

reference <- unique(fread(
  file.path(project_dir, "data", "era5_consistent",
            "era5_201506_consistent.csv"),
  select = c("site_id", "lon", "lat")
))[order(site_id)]
n_sites <- nrow(reference)
operators <- make_graph_operators(reference, h_factors)
coords <- project_coordinates_km(reference$lon, reference$lat)
distance <- as.matrix(dist(coords))
centering <- diag(n_sites) - matrix(1 / n_sites, n_sites, n_sites)

covariance_matrix <- function(range_km, nugget = 0.05) {
  (1 - nugget) * exp(-distance / range_km) + nugget * diag(n_sites)
}
sigma_short <- covariance_matrix(450)
sigma_weak <- covariance_matrix(525)
sigma_long <- covariance_matrix(650)
sigma_very_long <- covariance_matrix(800)
chol_short <- t(chol(sigma_short + 1e-8 * diag(n_sites)))
chol_weak <- t(chol(sigma_weak + 1e-8 * diag(n_sites)))
chol_long <- t(chol(sigma_long + 1e-8 * diag(n_sites)))
chol_very_long <- t(chol(sigma_very_long + 1e-8 * diag(n_sites)))
gradient <- as.numeric(scale(0.8 * coords[, 1] + coords[, 2]))

make_pair_operator <- function(W) {
  diag(W) <- 0
  W <- (W + t(W)) > 0
  storage.mode(W) <- "double"
  L <- diag(rowSums(W)) - W
  list(W = W, L = L, weight_sum = sum(W[upper.tri(W)]))
}

nearest_W <- matrix(0, n_sites, n_sites)
for (i in seq_len(n_sites)) {
  nearest <- order(distance[i, ])[2:5]
  nearest_W[i, nearest] <- 1
}
nearest_operator <- make_pair_operator(nearest_W)
binned_operator <- make_pair_operator(distance >= 400 & distance < 600)
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
variogram_names <- sprintf("variogram_bin_%d", seq_len(5))
moran_operator <- operators[[3]]

all_quadratic <- c(
  list(spatial_variance = list(
    L = centering,
    divisor = n_sites
  )),
  list(nearest_semivariance = list(
    L = nearest_operator$L,
    divisor = 2 * nearest_operator$weight_sum
  )),
  list(binned_semivariance = list(
    L = binned_operator$L,
    divisor = 2 * binned_operator$weight_sum
  )),
  setNames(lapply(variogram_operators, function(op) list(
    L = op$L, divisor = 2 * op$weight_sum
  )), variogram_names),
  setNames(lapply(operators, function(op) list(
    L = op$L, divisor = 2 * op$weight_sum
  )), graph_names)
)
quadratic_names <- names(all_quadratic)

expected_components <- function(sigma) {
  vapply(all_quadratic, function(op) {
    sum(diag(op$L %*% sigma)) / op$divisor
  }, numeric(1))
}
gradient_components <- vapply(all_quadratic, function(op) {
  as.numeric(crossprod(gradient, op$L %*% gradient) / op$divisor)
}, numeric(1))
expected_short <- expected_components(sigma_short)
expected_weak <- expected_components(sigma_weak)
expected_long <- expected_components(sigma_long)
expected_very_long <- expected_components(sigma_very_long)

scenario_table <- rbindlist(list(
  data.table(
    scenario = c("cyclic_base", "cyclic_3_records", "cyclic_12_records",
                 "cyclic_20_days", "cyclic_60_days", "cyclic_rho0",
                 "cyclic_rho075", "cyclic_t3", "linear_base",
                 "linear_rho075"),
    family = "null", mechanism = "null", strength = 0,
    null_type = c(rep("cyclic", 8), rep("linear", 2)),
    records = c(6, 3, 12, 6, 6, 6, 6, 6, 6, 6),
    days = c(30, 30, 30, 20, 60, 30, 30, 30, 30, 30),
    rho_field = c(0.45, 0.45, 0.45, 0.45, 0.45, 0, 0.75, 0.45, 0.45, 0.75),
    spatial_distribution = c(rep("gaussian", 7), "t3", "gaussian", "gaussian")
  ),
  data.table(
    scenario = paste0("amplitude_", c("weak", "weak_plus", "medium",
                                      "medium_plus", "strong")),
    family = "alternative", mechanism = "amplitude",
    strength = c(0.97, 0.955, 0.94, 0.92, 0.90), null_type = "linear",
    records = 6, days = 30, rho_field = 0.45,
    spatial_distribution = "gaussian"
  ),
  data.table(
    scenario = paste0("range_", c("weak", "weak_plus", "medium",
                                  "medium_plus", "strong")),
    family = "alternative", mechanism = "range",
    strength = c(525, 587, 650, 725, 800), null_type = "linear",
    records = 6, days = 30, rho_field = 0.45,
    spatial_distribution = "gaussian"
  ),
  data.table(
    scenario = paste0("gradient_", c("weak", "weak_plus", "medium",
                                     "medium_plus", "strong")),
    family = "alternative", mechanism = "gradient",
    strength = c(1.6, 1.5, 1.4, 1.3, 1.2), null_type = "linear",
    records = 6, days = 30, rho_field = 0.45,
    spatial_distribution = "gaussian"
  )
), fill = TRUE)

simulate_ar1 <- function(n, rho) {
  suppressWarnings(as.numeric(arima.sim(
    model = list(ar = rho), n = n, sd = sqrt(1 - rho^2)
  )))
}

scenario_chol <- function(row, extreme) {
  if (row$mechanism != "range" || !extreme) return(chol_short)
  if (row$strength < 600) return(chol_weak)
  if (row$strength < 750) chol_long else chol_very_long
}

simulate_record <- function(n, row) {
  mu <- simulate_ar1(n, 0.65)
  q75 <- quantile(mu, 0.75, names = FALSE)
  is_extreme <- mu >= q75
  field <- matrix(0, n_sites, n)

  if (row$null_type == "cyclic") {
    index_distance <- abs(outer(seq_len(n), seq_len(n), "-"))
    circular_distance <- pmin(index_distance, n - index_distance)
    temporal_covariance <- row$rho_field^circular_distance
    temporal_chol <- t(chol(temporal_covariance + 1e-8 * diag(n)))
    field <- chol_short %*% matrix(rnorm(n_sites * n), n_sites, n) %*%
      t(temporal_chol)
    if (row$spatial_distribution == "t3") {
      field <- sweep(field, 2, sqrt(3 / rchisq(n, df = 3)), "*")
    }
    field <- sweep(field, 2, colMeans(field), "-")
    field <- sweep(field, 2, mu, "+")
    return(list(mean = mu, field = field, extreme = is_extreme))
  }

  previous <- chol_short %*% rnorm(n_sites)
  for (tt in seq_len(n)) {
    innovation <- scenario_chol(row, is_extreme[tt]) %*% rnorm(n_sites)
    state <- row$rho_field * previous +
      sqrt(1 - row$rho_field^2) * innovation
    previous <- state
    anomaly <- state - mean(state)
    if (row$mechanism == "amplitude" && is_extreme[tt]) {
      anomaly <- row$strength * anomaly
    }
    if (row$mechanism == "gradient") {
      amplitude <- if (is_extreme[tt]) row$strength else 1.8
      anomaly <- anomaly + amplitude * gradient
    }
    field[, tt] <- mu[tt] + anomaly
  }
  list(mean = mu, field = field, extreme = is_extreme)
}

field_metrics <- function(field) {
  centered <- sweep(field, 2, colMeans(field), "-")
  out <- lapply(all_quadratic, function(op) {
    colSums(field * (op$L %*% field)) / op$divisor
  })
  names(out) <- quadratic_names
  variance <- out$spatial_variance
  W <- moran_operator$W
  W_total <- sum(W)
  cross <- colSums(centered * (W %*% centered))
  sumsquares <- colSums(centered^2)
  out$moran_i <- n_sites / W_total * cross / sumsquares
  graph_energy <- colSums(field * (moran_operator$L %*% field))
  out$geary_c <- (n_sites - 1) * graph_energy /
    (2 * moran_operator$weight_sum * sumsquares)
  as.data.table(out)
}

conditional_target <- function(record_objects, row) {
  record_targets <- vapply(record_objects, function(record) {
    ext_sum <- mod_sum <- setNames(
      numeric(length(quadratic_names)), quadratic_names)
    ext_n <- mod_n <- 0L
    expected_state <- expected_short
    for (tt in seq_along(record$extreme)) {
      extreme <- record$extreme[tt]
      innovation_expected <- expected_short
      if (row$mechanism == "range" && extreme) {
        innovation_expected <- if (row$strength < 600) expected_weak else
          if (row$strength < 750) expected_long else expected_very_long
      }
      expected_state <- row$rho_field^2 * expected_state +
        (1 - row$rho_field^2) * innovation_expected
      expected_output <- expected_state
      if (row$mechanism == "amplitude" && extreme) {
        expected_output <- row$strength^2 * expected_output
      }
      if (row$mechanism == "gradient") {
        amplitude <- if (extreme) row$strength else 1.8
        expected_output <- expected_output + amplitude^2 * gradient_components
      }
      if (extreme) {
        ext_sum <- ext_sum + expected_output
        ext_n <- ext_n + 1L
      } else {
        # Match the analysis: only the interquartile moderate days enter.
        mu <- record$mean
        q25 <- quantile(mu, 0.25, names = FALSE)
        if (mu[tt] >= q25) {
          mod_sum <- mod_sum + expected_output
          mod_n <- mod_n + 1L
        }
      }
    }
    ext_sum / ext_n / (mod_sum / mod_n) - 1
  }, numeric(length(quadratic_names)))
  relative <- rowMeans(record_targets)
  mean(relative[graph_names])
}

simulate_dataset <- function(row) {
  records <- lapply(seq_len(row$records), function(r) {
    simulate_record(row$days + as.integer(r %% 3L == 2L), row)
  })
  dat <- rbindlist(lapply(seq_along(records), function(r) {
    object <- records[[r]]
    cbind(
      data.table(file = sprintf("record_%02d", r),
                 date = seq_along(object$mean), wbt_mean = object$mean),
      field_metrics(object$field)
    )
  }))
  list(data = dat, target = conditional_target(records, row))
}

joint_shift_test <- function(dat, B, seed) {
  classified <- classify_mean_quantiles(dat)
  metric_names <- c(quadratic_names, "moran_i", "geary_c")
  positive_names <- quadratic_names
  strata <- lapply(unique(classified$file), function(f) {
    x <- classified[file == f][order(date)]
    list(regime = x$regime, outcomes = as.matrix(x[, ..metric_names]), n = nrow(x))
  })
  record_effect <- function(outcomes, regime) {
    ext_mean <- colMeans(outcomes[regime == "extreme", , drop = FALSE])
    mod_mean <- colMeans(outcomes[regime == "moderate", , drop = FALSE])
    effect <- ext_mean - mod_mean
    effect[positive_names] <- ext_mean[positive_names] /
      mod_mean[positive_names] - 1
    effect
  }
  observed <- rowMeans(vapply(strata, function(x) {
    record_effect(x$outcomes, x$regime)
  }, numeric(length(metric_names))))
  set.seed(seed)
  null <- matrix(0, B, length(metric_names),
                 dimnames = list(NULL, metric_names))
  for (x in strata) {
    is_ext <- x$regime == "extreme"
    is_mod <- x$regime == "moderate"
    effect_by_offset <- matrix(0, x$n, length(metric_names),
                               dimnames = list(NULL, metric_names))
    for (offset in 0:(x$n - 1L)) {
      index <- ((seq_len(x$n) - 1L + offset) %% x$n) + 1L
      shifted <- x$outcomes[index, , drop = FALSE]
      ext_mean <- colMeans(shifted[is_ext, , drop = FALSE])
      mod_mean <- colMeans(shifted[is_mod, , drop = FALSE])
      effect_by_offset[offset + 1L, ] <- ext_mean - mod_mean
      effect_by_offset[offset + 1L, positive_names] <-
        ext_mean[positive_names] / mod_mean[positive_names] - 1
    }
    sampled <- sample.int(x$n, B, replace = TRUE)
    null <- null + effect_by_offset[sampled, , drop = FALSE]
  }
  null <- null / length(strata)
  profile_observed <- mean(observed[graph_names])
  profile_null <- rowMeans(null[, graph_names, drop = FALSE])
  variogram_observed <- mean(observed[variogram_names])
  variogram_null <- rowMeans(null[, variogram_names, drop = FALSE])
  lower_p <- function(metric) (1 + sum(null[, metric] <= observed[metric])) / (B + 1)
  two_sided_p <- function(metric) {
    (1 + sum(abs(null[, metric]) >= abs(observed[metric]))) / (B + 1)
  }
  list(
    estimate_profile = profile_observed,
    p_graph = (1 + sum(profile_null <= profile_observed)) / (B + 1),
    p_variogram_profile =
      (1 + sum(variogram_null <= variogram_observed)) / (B + 1),
    p_variance = lower_p("spatial_variance"),
    p_nearest = lower_p("nearest_semivariance"),
    p_binned = lower_p("binned_semivariance"),
    p_moran = two_sided_p("moran_i"),
    p_geary = two_sided_p("geary_c"),
    effect_variance = observed["spatial_variance"],
    effect_nearest = observed["nearest_semivariance"],
    effect_binned = observed["binned_semivariance"],
    effect_local = observed[graph_names[1]],
    effect_broad = observed[graph_names[length(graph_names)]]
  )
}

set.seed(20260872L)
results <- vector("list", nrow(scenario_table) * n_sim)
counter <- 0L
for (j in seq_len(nrow(scenario_table))) {
  row <- scenario_table[j]
  cat(sprintf("Scenario %s (%d simulations)\n", row$scenario, n_sim))
  for (s in seq_len(n_sim)) {
    counter <- counter + 1L
    simulated <- simulate_dataset(row)
    test <- joint_shift_test(simulated$data, B, 20260872L + counter)
    results[[counter]] <- cbind(
      row,
      data.table(
        simulation = s,
        target_profile = simulated$target,
        estimate_profile = test$estimate_profile,
        estimation_error = test$estimate_profile - simulated$target,
        p_graph = test$p_graph,
        p_variogram_profile = test$p_variogram_profile,
        p_variance = test$p_variance,
        p_nearest = test$p_nearest,
        p_binned = test$p_binned,
        p_moran = test$p_moran,
        p_geary = test$p_geary,
        effect_variance = test$effect_variance,
        effect_nearest = test$effect_nearest,
        effect_binned = test$effect_binned,
        effect_local = test$effect_local,
        effect_broad = test$effect_broad
      )
    )
  }
}

results <- rbindlist(results, fill = TRUE)
summary <- results[, .(
  simulations = .N,
  target_profile = mean(target_profile),
  bias_profile = mean(estimation_error),
  rmse_profile = sqrt(mean(estimation_error^2)),
  rejection_graph = mean(p_graph <= 0.05),
  rejection_variogram_profile = mean(p_variogram_profile <= 0.05),
  rejection_variance = mean(p_variance <= 0.05),
  rejection_nearest = mean(p_nearest <= 0.05),
  rejection_binned = mean(p_binned <= 0.05),
  rejection_moran = mean(p_moran <= 0.05),
  rejection_geary = mean(p_geary <= 0.05),
  median_relative_variance = median(effect_variance),
  median_relative_nearest = median(effect_nearest),
  median_relative_binned = median(effect_binned),
  median_relative_local = median(effect_local),
  median_relative_broad = median(effect_broad)
), by = .(scenario, family, mechanism, strength, null_type, records, days,
          rho_field, spatial_distribution)]

fwrite(results, file.path(output_dir, "graph_simulation_results.csv"))
fwrite(summary, file.path(output_dir, "graph_simulation_summary.csv"))
print(summary)
