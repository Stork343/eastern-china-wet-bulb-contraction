###############################################################################
# Sensitivity analyses using the calibrated circular-shift test.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))

data_dir <- file.path(project_dir, "data", "era5_consistent")
output_dir <- file.path(project_dir, "output_corrected")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hourly <- load_consistent_era5(data_dir)
daily <- aggregate_daily_peak_fields(hourly)
primary_metrics <- compute_field_metrics(daily, h_factor = 0.5)

run_one <- function(metrics, metric, lower, upper, seed) {
  fit <- circular_shift_test(
    metrics, metric, B = 4999L, seed = seed, lower = lower, upper = upper
  )
  record <- record_level_inference(metrics, metric, lower = lower, upper = upper)
  data.table(
    metric = metric,
    lower_quantile = lower,
    upper_quantile = upper,
    delta = fit$observed,
    shift_p_value = fit$p_value,
    record_ci_lower = record$ci[1],
    record_ci_upper = record$ci[2]
  )
}

cat("Extreme-threshold sensitivity...\n")
threshold_results <- rbindlist(lapply(seq_along(c(0.70, 0.75, 0.80)), function(i) {
  upper <- c(0.70, 0.75, 0.80)[i]
  run_one(primary_metrics, "spatial_variance", 0.25, upper,
          20260820L + i)[, analysis := "extreme_threshold"]
}))

cat("Daily-field definition sensitivity...\n")
composite_metrics <- compute_field_metrics(
  aggregate_daily_composite_max(copy(hourly)), h_factor = 0.5
)
field_result <- run_one(
  composite_metrics, "spatial_variance", 0.25, 0.75, 20260829L
)[, `:=`(analysis = "field_definition",
         field_definition = "sitewise_daily_max_composite")]

cat("Wet-bulb algorithm sensitivity...\n")
stull_metrics <- compute_field_metrics(
  aggregate_daily_peak_fields(copy(hourly), wbt_column = "wbt_stull"),
  h_factor = 0.5
)
wbt_method_result <- run_one(
  stull_metrics, "spatial_variance", 0.25, 0.75, 20260829L
)[, `:=`(analysis = "wbt_algorithm", wbt_algorithm = "Stull_2011")]

cat("Kernel-bandwidth sensitivity...\n")
bandwidth_results <- rbindlist(lapply(seq_along(c(0.25, 0.5, 1, 2)), function(i) {
  factor <- c(0.25, 0.5, 1, 2)[i]
  metrics <- compute_field_metrics(daily, h_factor = factor)
  result <- run_one(metrics, "lambda_weighted", 0.25, 0.75,
                    20260830L + i)
  result[, `:=`(
    analysis = "kernel_bandwidth",
    h_factor = factor,
    bandwidth_km = attr(metrics, "bandwidth_km")
  )]
  result
}))

all_results <- rbindlist(
  list(threshold_results, field_result, wbt_method_result, bandwidth_results),
  use.names = TRUE, fill = TRUE
)
fwrite(all_results, file.path(output_dir, "corrected_sensitivity.csv"))
print(all_results, digits = 4)
