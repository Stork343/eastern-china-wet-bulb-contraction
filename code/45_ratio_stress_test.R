###############################################################################
# Heavy-tail stress test for raw, log, and bounded ratio effect measures.
#
# The design is frozen in EXTENSION_ANALYSIS_PROTOCOL.md.  It reuses the
# cyclic-t3 day-level field mechanism from 18_graph_esh_simulation.R.  Each
# replication generates one six-record heavy-tailed field panel.  The null and
# -7% alternatives are paired versions of that panel, and all three effect
# measures use exactly the same product cyclic shifts.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))

output_dir <- file.path(project_dir, "output_extension_methods")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_sim <- as.integer(Sys.getenv("ESH_RATIO_STRESS_N", unset = "2000"))
B <- as.integer(Sys.getenv("ESH_RATIO_STRESS_B", unset = "999"))
seed <- 20260811L
h_factors <- c(0.125, 0.25, 0.5, 1, 2)
record_lengths <- c(30L, 31L, 30L, 30L, 31L, 30L)
rho_mean <- 0.65
rho_field <- 0.45
t_df <- 3
alternative_ratio <- 0.93
measure_names <- c("raw_ratio", "log_ratio", "bounded_symmetric")

reference_path <- file.path(
  project_dir, "data", "era5_consistent", "era5_201506_consistent.csv"
)
reference <- unique(fread(reference_path, select = c("site_id", "lon", "lat")))[
  order(site_id)
]
n_sites <- nrow(reference)
if (n_sites != 121L) stop("Expected 121 reference sites")

coords <- project_coordinates_km(reference$lon, reference$lat)
distance <- as.matrix(dist(coords))
sigma <- (1 - 0.05) * exp(-distance / 450) + 0.05 * diag(n_sites)
spatial_chol <- t(chol(sigma + 1e-8 * diag(n_sites)))
operators <- make_graph_operators(reference, h_factors)
bandwidths <- vapply(operators, `[[`, numeric(1), "bandwidth_km")

# Applying C'LC to standard-normal spatial coordinates removes one matrix
# multiplication per record while producing the same quadratic forms as the
# original spatial field generator.
standard_operators <- lapply(operators, function(operator) {
  list(
    A = t(spatial_chol) %*% operator$L %*% spatial_chol,
    divisor = 2 * operator$weight_sum
  )
})

simulate_ar1 <- function(n, rho) {
  suppressWarnings(as.numeric(arima.sim(
    model = list(ar = rho), n = n, sd = sqrt(1 - rho^2)
  )))
}

classify_regime <- function(mu) {
  thresholds <- quantile(mu, c(0.25, 0.75), names = FALSE, type = 7)
  fifelse(mu >= thresholds[2], "high",
          fifelse(mu >= thresholds[1], "middle", "low"))
}

temporal_cholesky <- lapply(sort(unique(record_lengths)), function(n) {
  index_distance <- abs(outer(seq_len(n), seq_len(n), "-"))
  circular_distance <- pmin(index_distance, n - index_distance)
  covariance <- rho_field^circular_distance
  t(chol(covariance + 1e-8 * diag(n)))
})
names(temporal_cholesky) <- as.character(sort(unique(record_lengths)))

simulate_base_record <- function(n) {
  mu <- simulate_ar1(n, rho_mean)
  regime <- classify_regime(mu)
  temporal_chol <- temporal_cholesky[[as.character(n)]]
  standard_field <- matrix(rnorm(n_sites * n), n_sites, n) %*%
    t(temporal_chol)
  # This common day multiplier exactly matches the cyclic_t3 mechanism in
  # 18_graph_esh_simulation.R.  Its square makes graph energy heavy tailed.
  day_multiplier <- sqrt(t_df / rchisq(n, df = t_df))
  metrics <- vapply(standard_operators, function(operator) {
    colSums(standard_field * (operator$A %*% standard_field)) /
      operator$divisor * day_multiplier^2
  }, numeric(n))
  colnames(metrics) <- sprintf("graph_%d", seq_along(h_factors))
  list(regime = regime, metrics = metrics, n = n)
}

transform_ratio <- function(ratio, measure) {
  if (measure == "raw_ratio") return(ratio - 1)
  if (measure == "log_ratio") return(log(ratio))
  if (measure == "bounded_symmetric") return(2 * (ratio - 1) / (ratio + 1))
  stop("Unknown effect measure: ", measure)
}

measure_target <- function(ratio, measure) transform_ratio(ratio, measure)

effect_by_offset <- function(metrics, regime, measure) {
  n <- nrow(metrics)
  high <- regime == "high"
  middle <- regime == "middle"
  result <- matrix(NA_real_, n, ncol(metrics))
  for (offset in 0:(n - 1L)) {
    index <- ((seq_len(n) - 1L + offset) %% n) + 1L
    shifted <- metrics[index, , drop = FALSE]
    ratio <- colMeans(shifted[high, , drop = FALSE]) /
      colMeans(shifted[middle, , drop = FALSE])
    result[offset + 1L, ] <- transform_ratio(ratio, measure)
  }
  result
}

paired_tests <- function(base_records) {
  # Draw one product-shift matrix and reuse it for null/alternative, all five
  # scales, and all three effect measures.
  sampled_offsets <- lapply(base_records, function(record) {
    sample.int(record$n, B, replace = TRUE)
  })
  scenario_records <- list(
    null = base_records,
    alternative_minus_7pct = lapply(base_records, function(record) {
      copy_record <- record
      copy_record$metrics <- record$metrics
      copy_record$metrics[record$regime == "high", ] <-
        alternative_ratio * copy_record$metrics[record$regime == "high", ]
      copy_record
    })
  )

  rows <- list()
  counter <- 0L
  for (scenario in names(scenario_records)) {
    records <- scenario_records[[scenario]]
    true_ratio <- if (scenario == "null") 1 else alternative_ratio
    for (measure in measure_names) {
      lookup <- lapply(records, function(record) {
        effect_by_offset(record$metrics, record$regime, measure)
      })
      observed_scale <- Reduce(`+`, lapply(lookup, function(x) x[1, ])) /
        length(lookup)
      null_scale <- matrix(0, B, length(h_factors))
      for (r in seq_along(lookup)) {
        null_scale <- null_scale + lookup[[r]][sampled_offsets[[r]], , drop = FALSE]
      }
      null_scale <- null_scale / length(lookup)
      observed_profile <- mean(observed_scale)
      null_profile <- rowMeans(null_scale)
      p_value <- (1 + sum(null_profile <= observed_profile)) / (B + 1)
      counter <- counter + 1L
      rows[[counter]] <- data.table(
        scenario = scenario,
        measure = measure,
        target = measure_target(true_ratio, measure),
        estimate = observed_profile,
        p_value = p_value,
        reject_05 = p_value <= 0.05
      )
    }
  }
  rbindlist(rows)
}

set.seed(seed)
result_rows <- vector("list", n_sim)
cat(sprintf(
  "Running %d paired heavy-tail replications with %d product shifts each...\n",
  n_sim, B
))
for (simulation in seq_len(n_sim)) {
  records <- lapply(record_lengths, simulate_base_record)
  result_rows[[simulation]] <- paired_tests(records)[, simulation := simulation]
  if (simulation %% 50L == 0L || simulation == n_sim) {
    cat(sprintf("  completed %d/%d\n", simulation, n_sim))
  }
}
results <- rbindlist(result_rows)
setcolorder(results, c("simulation", "scenario", "measure", "target",
                       "estimate", "p_value", "reject_05"))

summary <- results[, .(
  replications = .N,
  target = unique(target),
  mean_estimate = mean(estimate),
  bias = mean(estimate - target),
  rmse = sqrt(mean((estimate - target)^2)),
  rejection_rate_05 = mean(reject_05),
  monte_carlo_se_rejection = sqrt(
    mean(reject_05) * (1 - mean(reject_05)) / .N
  ),
  median_p_value = median(p_value)
), by = .(scenario, measure)]

comparison <- dcast(
  summary,
  measure ~ scenario,
  value.var = c("target", "mean_estimate", "bias", "rmse",
                "rejection_rate_05", "monte_carlo_se_rejection")
)

results_path <- file.path(output_dir, "extension_ratio_stress_results.csv.gz")
summary_path <- file.path(output_dir, "extension_ratio_stress_summary.csv")
comparison_path <- file.path(output_dir, "extension_ratio_stress_comparison.csv")
fwrite(results, results_path)
fwrite(summary, summary_path)
fwrite(comparison, comparison_path)

label <- c(
  raw_ratio = "Raw ratio",
  log_ratio = "Log ratio",
  bounded_symmetric = "Bounded symmetric"
)
tex_path <- file.path(output_dir, "extension_ratio_stress_table.tex")
tex <- c(
  "% Generated by code/45_ratio_stress_test.R",
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  paste0("Measure & Null bias & Null RMSE & Size & ",
         "$-7\\%$ bias & $-7\\%$ RMSE & Power \\\\"),
  "\\midrule"
)
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i]
  tex <- c(tex, sprintf(
    "%s & %.4f & %.4f & %.3f & %.4f & %.4f & %.3f \\\\",
    label[[row$measure]], row$bias_null, row$rmse_null,
    row$rejection_rate_05_null, row$bias_alternative_minus_7pct,
    row$rmse_alternative_minus_7pct,
    row$rejection_rate_05_alternative_minus_7pct
  ))
}
tex <- c(
  tex, "\\bottomrule", "\\end{tabular}",
  paste0("% Heavy-tailed cyclic-t3 DGP; ", n_sim,
         " paired replications per scenario; B=", B,
         " shared product shifts; seed=", seed, ".")
)
writeLines(tex, tex_path, useBytes = TRUE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
outputs <- c(results_path, summary_path, comparison_path, tex_path)
audit <- list(
  protocol = "EXTENSION_ANALYSIS_PROTOCOL.md",
  script = file.path("code", basename(script_path)),
  analysis_role = "post-analysis ratio stress test",
  seed = seed,
  replications_per_cell = n_sim,
  product_shift_draws = B,
  scenarios = c("null", "alternative_minus_7pct"),
  effect_measures = measure_names,
  paired_scenarios = TRUE,
  shared_product_shifts_across_scenarios_scales_measures = TRUE,
  site_count = n_sites,
  record_lengths = record_lengths,
  fixed_bandwidths_km = bandwidths,
  dgp = list(
    spatial_covariance = "0.95*exp(-distance_km/450)+0.05*I",
    circular_temporal_correlation = rho_field,
    regional_mean_ar1 = rho_mean,
    common_day_scale = "sqrt(3/chi_square_3)",
    alternative_high_day_graph_energy_multiplier = alternative_ratio
  ),
  inputs = list(list(
    path = file.path("data", "era5_consistent", basename(reference_path)),
    bytes = unname(file.info(reference_path)$size),
    sha256 = sha256_file(reference_path)
  )),
  outputs = lapply(outputs, function(path) list(
    path = file.path("output_extension_methods", basename(path)),
    bytes = unname(file.info(path)$size),
    sha256 = sha256_file(path)
  ))
)
jsonlite::write_json(
  audit,
  file.path(output_dir, "extension_ratio_stress_audit.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)

print(summary)
