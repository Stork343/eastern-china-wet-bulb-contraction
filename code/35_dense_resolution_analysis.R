###############################################################################
# Spatial-resolution sensitivity on the nested 465-site ERA5-Land grid.
#
# Configurations:
#   1. Reproduce the primary analysis on the embedded original 121 sites.
#   2. Use all 465 sites while freezing primary peak times and day labels.
#   3. Use all 465 sites and recompute peak times and labels on the dense grid.
#
# Absolute graph bandwidths remain fixed at the primary 126--2013 km values.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))

input_dir <- file.path(project_dir, "data", "era5_dense", "daily_fields")
output_dir <- file.path(project_dir, "output_dense")
manifest_file <- file.path(project_dir, "data", "grid",
                           "eastern_china_dense_sites.csv")
primary_output <- file.path(project_dir, "output_confirmatory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = "^era5_land_[0-9]{4}_jja_dense_daily_fields\\.csv\\.gz$",
  full.names = TRUE
)
if (length(files) != 35L) stop("Expected 35 dense daily-field files")

keep <- c("analysis_role", "year", "month", "record_id", "analysis_date",
          "analysis_definition", "site_id", "requested_lon", "requested_lat",
          "label_mean_wbt", "regional_mean_wbt", "wbt")
fields <- rbindlist(lapply(files, function(path) fread(path, select = keep)))
fields[, `:=`(
  file = as.character(record_id),
  date = as.IDate(analysis_date),
  lon = requested_lon,
  lat = requested_lat
)]
setorder(fields, analysis_definition, file, date, site_id)
if (fields[, uniqueN(site_id)] != 465L) stop("Expected 465 dense sites")
counts <- fields[, .N, by = .(analysis_definition, file, date)]
if (counts[, any(N != 465L)]) stop("Incomplete dense daily field")

manifest <- fread(manifest_file)
if (nrow(manifest) != 465L || manifest[, sum(is_original_site)] != 121L) {
  stop("Invalid dense manifest")
}

primary_metadata <- fread(file.path(
  primary_output, "confirmatory_graph_metadata.csv"
))[definition_index == 1L]
primary_metadata <- unique(primary_metadata[, .(
  metric, h_factor, bandwidth_km
)])
setorder(primary_metadata, h_factor)
if (nrow(primary_metadata) != 5L) stop("Expected five primary bandwidths")

make_fixed_operators <- function(sites, metadata) {
  sites <- copy(sites)[order(site_id)]
  coords <- project_coordinates_km(sites$lon, sites$lat)
  distance <- as.matrix(dist(coords))
  lapply(seq_len(nrow(metadata)), function(k) {
    h <- metadata$bandwidth_km[k]
    W <- exp(-(distance^2) / (2 * h^2))
    diag(W) <- 0
    L <- diag(rowSums(W)) - W
    list(
      metric = metadata$metric[k],
      h_factor = metadata$h_factor[k],
      bandwidth_km = h,
      W = W,
      L = L,
      weight_sum = sum(W[upper.tri(W)])
    )
  })
}

compute_configuration <- function(definition, site_ids, configuration) {
  x <- fields[
    analysis_definition == definition & site_id %chin% site_ids
  ]
  sites <- unique(x[, .(site_id, lon, lat)])[order(site_id)]
  if (nrow(sites) != length(site_ids)) stop("Configuration site mismatch")
  operators <- make_fixed_operators(sites, primary_metadata)

  metadata <- unique(x[, .(
    file, date, year, month, analysis_role, label_mean_wbt,
    regional_mean_wbt
  )])
  if (metadata[, anyDuplicated(paste(file, date))]) {
    stop("Nonunique daily metadata in ", configuration)
  }
  wide <- dcast(x[, .(file, date, site_id, wbt)],
                file + date ~ site_id, value.var = "wbt")
  value_columns <- as.character(sites$site_id)
  if (length(setdiff(value_columns, names(wide)))) {
    stop("Wide field matrix is missing sites")
  }
  Y <- as.matrix(wide[, ..value_columns])
  metrics <- wide[, .(file, date)]
  for (op in operators) {
    metrics[, (op$metric) := rowSums((Y %*% op$L) * Y) /
      (2 * op$weight_sum)]
  }
  metrics <- merge(metrics, metadata, by = c("file", "date"), all.x = TRUE)
  metrics[, wbt_mean := label_mean_wbt]
  classified <- classify_mean_quantiles(metrics, lower = 0.25, upper = 0.75)

  metric_names <- primary_metadata$metric
  record_effects <- rbindlist(lapply(metric_names, function(metric) {
    classified[, .(
      mean_extreme = mean(get(metric)[regime == "extreme"]),
      mean_moderate = mean(get(metric)[regime == "moderate"]),
      relative_effect = mean(get(metric)[regime == "extreme"]) /
        mean(get(metric)[regime == "moderate"]) - 1,
      absolute_effect = mean(get(metric)[regime == "extreme"]) -
        mean(get(metric)[regime == "moderate"]),
      n_extreme = sum(regime == "extreme"),
      n_moderate = sum(regime == "moderate")
    ), by = .(file, year, month, analysis_role)][, metric := metric]
  }))
  record_effects[, configuration := configuration]

  year_effects <- record_effects[, .(
    yearly_relative_effect = mean(relative_effect),
    months = .N
  ), by = .(configuration, analysis_role, year, metric)]
  if (year_effects[, any(months != 3L)]) stop("Incomplete dense yearly effects")
  profile_years <- year_effects[, .(
    yearly_profile_effect = mean(yearly_relative_effect),
    scales = .N
  ), by = .(configuration, analysis_role, year)]
  if (profile_years[, any(scales != 5L)]) stop("Incomplete dense profiles")

  list(
    metrics = classified,
    records = record_effects,
    years = year_effects,
    profiles = profile_years,
    operators = operators
  )
}

original_dense_ids <- manifest[is_original_site == TRUE, dense_site_id]
fits <- list(
  primary_121_reproduction = compute_configuration(
    "primary_grid_peak", original_dense_ids, "primary_121_reproduction"
  ),
  dense_465_fixed_labels = compute_configuration(
    "primary_grid_peak", manifest$dense_site_id, "dense_465_fixed_labels"
  ),
  dense_465_recomputed = compute_configuration(
    "dense_grid_peak", manifest$dense_site_id, "dense_465_recomputed"
  )
)

record_effects <- rbindlist(lapply(fits, `[[`, "records"))
year_effects <- rbindlist(lapply(fits, `[[`, "years"))
profile_years <- rbindlist(lapply(fits, `[[`, "profiles"))

one_sample_summary <- function(values, alternative = "less", hac_lag = 2L) {
  n <- length(values)
  estimate <- mean(values)
  standard_error <- sd(values) / sqrt(n)
  statistic <- estimate / standard_error
  t_p_value <- if (alternative == "less") pt(statistic, df = n - 1L) else {
    2 * pt(-abs(statistic), df = n - 1L)
  }
  centered <- values - estimate
  long_run_variance <- sum(centered^2) / n
  for (ell in seq_len(hac_lag)) {
    covariance <- sum(centered[(ell + 1L):n] *
                        centered[seq_len(n - ell)]) / n
    long_run_variance <- long_run_variance +
      2 * (1 - ell / (hac_lag + 1)) * covariance
  }
  hac_standard_error <- sqrt(long_run_variance / n)
  hac_statistic <- estimate / hac_standard_error
  hac_p_value <- if (alternative == "less") {
    pt(hac_statistic, df = n - 1L)
  } else {
    2 * pt(-abs(hac_statistic), df = n - 1L)
  }
  negative_years <- sum(values < 0)
  sign_p_value <- pbinom(negative_years - 1L, size = n, prob = 0.5,
                         lower.tail = FALSE)
  critical <- qt(0.975, df = n - 1L)
  data.table(
    years = n,
    estimate = estimate,
    standard_error = standard_error,
    t_statistic = statistic,
    t_p_value = t_p_value,
    ci_lower = estimate - critical * standard_error,
    ci_upper = estimate + critical * standard_error,
    hac_lag = hac_lag,
    hac_standard_error = hac_standard_error,
    hac_statistic = hac_statistic,
    hac_p_value = hac_p_value,
    hac_ci_lower = estimate - critical * hac_standard_error,
    hac_ci_upper = estimate + critical * hac_standard_error,
    negative_years = negative_years,
    negative_fraction = negative_years / n,
    sign_p_value = sign_p_value,
    consistency = t_p_value <= 0.025 & hac_p_value <= 0.025 &
      sign_p_value <= 0.025
  )
}

primary_results <- profile_years[analysis_role == "confirmatory",
  one_sample_summary(yearly_profile_effect), by = configuration]
scale_results <- year_effects[analysis_role == "confirmatory",
  one_sample_summary(yearly_relative_effect), by = .(configuration, metric)]
scale_results[, `:=`(
  t_holm_p = p.adjust(t_p_value, method = "holm"),
  hac_holm_p = p.adjust(hac_p_value, method = "holm"),
  sign_holm_p = p.adjust(sign_p_value, method = "holm")
), by = configuration]
scale_results[, scale_consistency := t_holm_p <= 0.025 &
                hac_holm_p <= 0.025 & sign_holm_p <= 0.025]
scale_results <- merge(scale_results, primary_metadata, by = "metric")

expected_years <- fread(file.path(
  primary_output, "confirmatory_year_profile_effects.csv"
))[day_definition == "utc" & analysis_role == "confirmatory",
   .(year, expected_yearly_profile_effect = yearly_profile_effect)]
reproduced_years <- profile_years[
  configuration == "primary_121_reproduction" &
    analysis_role == "confirmatory",
  .(year, reproduced_yearly_profile_effect = yearly_profile_effect)
]
reproduction_check <- merge(expected_years, reproduced_years, by = "year")
reproduction_check[, absolute_error := abs(
  reproduced_yearly_profile_effect - expected_yearly_profile_effect
)]
if (nrow(reproduction_check) != 33L ||
    reproduction_check[, max(absolute_error)] > 1e-9) {
  stop("Embedded 121-site analysis does not reproduce the primary result")
}

primary_reference <- fread(file.path(
  primary_output, "confirmatory_primary_results.csv"
))[day_definition == "utc", estimate]
primary_results[, `:=`(
  primary_reference = primary_reference,
  difference_from_primary = estimate - primary_reference
)]

fwrite(record_effects, file.path(output_dir, "dense_record_effects.csv"))
fwrite(year_effects, file.path(output_dir, "dense_year_scale_effects.csv"))
fwrite(profile_years, file.path(output_dir, "dense_year_profile_effects.csv"))
fwrite(primary_results, file.path(output_dir, "dense_primary_results.csv"))
fwrite(scale_results, file.path(output_dir, "dense_scale_results.csv"))
fwrite(reproduction_check,
       file.path(output_dir, "dense_primary_reproduction_check.csv"))

cat("Dense-grid primary summaries:\n")
print(primary_results[, .(
  configuration, years, estimate, ci_lower, ci_upper, t_p_value,
  hac_p_value, negative_years, sign_p_value, consistency,
  difference_from_primary
)])
cat(sprintf("Maximum embedded-grid reproduction error: %.3e\n",
            reproduction_check[, max(absolute_error)]))
