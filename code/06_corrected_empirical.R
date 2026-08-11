###############################################################################
# Corrected empirical analysis for extreme spatial homogenization.
#
# Primary choices:
#   - one observation is the synchronous field at the hour of the daily
#     maximum regional-mean WBT;
#   - all months use the same 121-site grid;
#   - extreme labels use spatial mean WBT only, never spatial variance;
#   - thresholds are estimated within each month-year record;
#   - circular shifts within month-years preserve serial dependence under the
#     null while breaking mean--variance alignment;
#   - the primary outcome is spatial variance; weighted lambda is sensitivity.
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

cat("Loading corrected ERA5 data...\n")
hourly <- load_consistent_era5(data_dir)
daily <- aggregate_daily_peak_fields(hourly)
metrics <- compute_field_metrics(daily, h_factor = 0.5)

cat(sprintf("Daily fields: %d; fixed sites per field: %d\n",
            nrow(metrics), unique(metrics$n_sites)))
cat(sprintf("Weighted-kernel bandwidth: %.1f km\n",
            attr(metrics, "bandwidth_km")))

cat("Running within-record circular-shift tests...\n")
primary <- circular_shift_test(
  metrics, metric = "spatial_variance", B = 9999L, seed = 20260802L
)
weighted <- circular_shift_test(
  metrics, metric = "lambda_weighted", B = 9999L, seed = 20260803L
)
primary_record <- record_level_inference(metrics, "spatial_variance")
weighted_record <- record_level_inference(metrics, "lambda_weighted")
continuous <- circular_shift_slope_test(
  metrics, metric = "spatial_variance", B = 9999L, seed = 20260804L
)

summary <- data.table(
  metric = c("Spatial variance", "Weighted leading eigenvalue"),
  delta = c(primary$observed, weighted$observed),
  record_mean_delta = c(primary_record$estimate, weighted_record$estimate),
  record_ci_lower = c(primary_record$ci[1], weighted_record$ci[1]),
  record_ci_upper = c(primary_record$ci[2], weighted_record$ci[2]),
  record_t_p_value = c(primary_record$p_value_lower,
                       weighted_record$p_value_lower),
  shift_p_value = c(primary$p_value, weighted$p_value),
  n_records = 6L,
  shift_B = c(primary$B, weighted$B)
)

monthly_primary <- monthly_effects(metrics, "spatial_variance")
monthly_weighted <- monthly_effects(metrics, "lambda_weighted")
setnames(monthly_weighted,
         c("mean_extreme", "mean_moderate", "delta"),
         c("weighted_mean_extreme", "weighted_mean_moderate", "weighted_delta"))
monthly <- merge(monthly_primary, monthly_weighted[
  , .(file, weighted_mean_extreme, weighted_mean_moderate, weighted_delta)
], by = "file")

classified <- primary$classified
counts <- classified[, .N, by = regime][order(regime)]

fwrite(summary, file.path(output_dir, "corrected_esh_summary.csv"))
fwrite(monthly, file.path(output_dir, "corrected_monthly_effects.csv"))
fwrite(metrics, file.path(output_dir, "corrected_daily_metrics.csv"))
fwrite(data.table(
  model = "Spatial variance ~ within-record centered mean WBT + record FE",
  slope_per_degree_c = continuous$slope,
  shift_p_value = continuous$p_value,
  shift_B = continuous$B
), file.path(output_dir, "corrected_continuous_association.csv"))
saveRDS(list(primary = primary, weighted = weighted,
             primary_record = primary_record,
             weighted_record = weighted_record),
        file.path(output_dir, "corrected_inference_results.rds"))

cat("\nClassification counts:\n")
print(counts)
cat("\nCorrected ESH inference (shift test and record-level interval):\n")
print(summary, digits = 4)
cat("\nMonth-year effects:\n")
print(monthly, digits = 4)
cat("\nContinuous fixed-effect slope:\n")
print(data.table(slope_per_degree_c = continuous$slope,
                 shift_p_value = continuous$p_value), digits = 4)
cat("\nOutputs written to: ", output_dir, "\n", sep = "")
