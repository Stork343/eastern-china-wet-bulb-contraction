###############################################################################
# Multi-scale graph-dispersion analysis for the corrected ERA5-Land fields.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))

data_dir <- file.path(project_dir, "data", "era5_consistent")
output_dir <- file.path(project_dir, "output_corrected")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

h_factors <- c(0.125, 0.25, 0.5, 1, 2)
hourly <- load_consistent_era5(data_dir)
daily <- aggregate_daily_peak_fields(hourly)
base_metrics <- compute_field_metrics(daily, h_factor = 0.5)
profile <- compute_graph_profile(daily, h_factors)
metrics <- merge(base_metrics, profile$metrics, by = c("file", "date"))
metric_names <- profile$metadata$metric

fit <- circular_shift_profile_test(
  metrics, metric_names, B = 9999L, seed = 20260870L
)

results <- merge(profile$metadata, fit$effects, by = "metric")
record_results <- rbindlist(lapply(metric_names, function(metric) {
  by_record <- fit$classified[, .(
    relative_delta_record = mean(get(metric)[regime == "extreme"]) /
      mean(get(metric)[regime == "moderate"]) - 1
  ), by = file]
  estimate <- mean(by_record$relative_delta_record)
  se <- sd(by_record$relative_delta_record) / sqrt(nrow(by_record))
  critical <- qt(0.975, df = nrow(by_record) - 1L)
  data.table(
    metric = metric,
    record_mean_relative = estimate,
    record_ci_lower = estimate - critical * se,
    record_ci_upper = estimate + critical * se
  )
}))
results <- merge(results, record_results, by = "metric")
results[, `:=`(
  omnibus_statistic = fit$omnibus_statistic,
  omnibus_p = fit$omnibus_p,
  shift_B = fit$B
)]

fwrite(metrics, file.path(output_dir, "graph_daily_metrics.csv"))
fwrite(results, file.path(output_dir, "graph_profile_results.csv"))
saveRDS(fit, file.path(output_dir, "graph_profile_shift_test.rds"))

cat("Multi-scale graph-dispersion results:\n")
print(results, digits = 4)
cat(sprintf("\nDirected profile omnibus: mean relative effect = %.3f, p = %.4g\n",
            fit$omnibus_statistic, fit$omnibus_p))

