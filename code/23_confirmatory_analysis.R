###############################################################################
# Confirmatory multi-year graph-dispersion analysis.
#
# The year is the replication unit. Discovery years 2015 and 2022 are retained
# in descriptive output but excluded from confirmatory tests.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))

input_dir <- Sys.getenv(
  "ESH_CONFIRMATORY_INPUT_DIR",
  unset = file.path(project_dir, "data", "era5_confirmatory", "daily_fields")
)
output_dir <- Sys.getenv(
  "ESH_CONFIRMATORY_OUTPUT_DIR",
  unset = file.path(project_dir, "output_confirmatory")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir,
                    pattern = "^era5_land_[0-9]{4}_jja_daily_fields\\.csv\\.gz$",
                    full.names = TRUE)
if (!length(files)) stop("No confirmatory daily-field files in ", input_dir)

fields <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
required <- c("analysis_role", "year", "month", "record_id", "analysis_date",
              "day_definition", "site_id", "requested_lon", "requested_lat",
              "regional_mean_wbt", "wbt")
if (length(missing <- setdiff(required, names(fields)))) {
  stop("Missing daily-field columns: ", paste(missing, collapse = ", "))
}

fields[, `:=`(
  file = as.character(record_id),
  date = as.IDate(analysis_date),
  lon = requested_lon,
  lat = requested_lat
)]
if (fields[, anyNA(date)]) stop("Failed to parse analysis dates")

h_factors <- c(0.125, 0.25, 0.5, 1, 2)
metric_names <- sprintf("graph_h_%s", format(h_factors, trim = TRUE))
discovery_years <- c(2015L, 2022L)

one_definition <- function(definition, lower = 0.25, upper = 0.75) {
  x <- fields[day_definition == definition]
  years <- sort(unique(x$year))
  if (!length(years)) stop("No data for day definition: ", definition)

  counts <- x[, .(sites = uniqueN(site_id)), by = .(file, date)]
  if (counts[, any(sites != 121L)]) stop("Incomplete spatial fields in ", definition)
  day_counts <- unique(x[, .(year, date)])[, .N, by = year]
  if (day_counts[, any(N != 92L)]) stop("Expected 92 days per year in ", definition)

  daily <- x[, .(file, date, site_id, lon, lat, wbt)]
  profile <- compute_graph_profile(daily, h_factors)
  base <- compute_field_metrics(daily, h_factor = 0.5)
  metrics <- merge(base[, .(file, date, wbt_mean, spatial_variance)],
                   profile$metrics, by = c("file", "date"))
  classified <- classify_mean_quantiles(metrics, lower = lower, upper = upper)

  record_effects <- rbindlist(lapply(c("spatial_variance", metric_names),
                                     function(metric) {
    classified[, .(
      mean_extreme = mean(get(metric)[regime == "extreme"]),
      mean_moderate = mean(get(metric)[regime == "moderate"]),
      relative_effect = mean(get(metric)[regime == "extreme"]) /
        mean(get(metric)[regime == "moderate"]) - 1,
      absolute_effect = mean(get(metric)[regime == "extreme"]) -
        mean(get(metric)[regime == "moderate"]),
      n_extreme = sum(regime == "extreme"),
      n_moderate = sum(regime == "moderate")
    ), by = file][, metric := metric]
  }))
  record_effects[, `:=`(
    year = as.integer(substr(file, 1, 4)),
    month = as.integer(substr(file, 5, 6)),
    day_definition = definition
  )]
  record_effects[, analysis_role := fifelse(year %in% discovery_years,
                                             "discovery", "confirmatory")]

  year_effects <- record_effects[, .(
    yearly_relative_effect = mean(relative_effect),
    months = .N
  ), by = .(day_definition, analysis_role, year, metric)]
  if (year_effects[, any(months != 3L)]) stop("Incomplete yearly record effects")

  profile_years <- year_effects[metric %chin% metric_names, .(
    yearly_profile_effect = mean(yearly_relative_effect),
    scales = .N
  ), by = .(day_definition, analysis_role, year)]
  if (profile_years[, any(scales != length(metric_names))]) {
    stop("Incomplete yearly scale profiles")
  }

  metrics[, `:=`(
    profile_dispersion = rowMeans(.SD),
    year = as.integer(substr(file, 1, 4)),
    month = as.integer(substr(file, 5, 6))
  ), .SDcols = metric_names]
  metrics[, analysis_role := fifelse(year %in% discovery_years,
                                      "discovery", "confirmatory")]

  mechanism_names <- intersect(c("u10", "v10", "swvl1", "ssrd"), names(x))
  mechanisms <- NULL
  if (length(mechanism_names)) {
    mechanisms <- x[, c(
      list(wind_speed_mean = if (all(c("u10", "v10") %chin% names(x))) {
        mean(sqrt(u10^2 + v10^2))
      } else NA_real_),
      lapply(.SD, mean)
    ), by = .(file, date), .SDcols = intersect(c("swvl1", "ssrd"), names(x))]
  }

  list(record = record_effects, year = year_effects, profile = profile_years,
       metrics = metrics, mechanisms = mechanisms,
       metadata = profile$metadata, lower = lower, upper = upper)
}

definitions <- intersect(c("utc", "utc_plus_8"), unique(fields$day_definition))
fits <- lapply(definitions, one_definition)
record_effects <- rbindlist(lapply(fits, `[[`, "record"))
year_effects <- rbindlist(lapply(fits, `[[`, "year"))
profile_years <- rbindlist(lapply(fits, `[[`, "profile"))

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
    degrees_freedom = n - 1L,
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
    consistency_threshold = 0.025,
    confirmatory_consistency = t_p_value <= 0.025 &
      hac_p_value <= 0.025 & sign_p_value <= 0.025
  )
}

primary <- profile_years[analysis_role == "confirmatory", {
  summary <- one_sample_summary(yearly_profile_effect)
  summary[, `:=`(
    loo_min = min(vapply(seq_along(yearly_profile_effect), function(i) {
      mean(yearly_profile_effect[-i])
    }, numeric(1))),
    loo_max = max(vapply(seq_along(yearly_profile_effect), function(i) {
      mean(yearly_profile_effect[-i])
    }, numeric(1)))
  )]
  summary
}, by = day_definition]

scale_results <- year_effects[
  analysis_role == "confirmatory" & metric %chin% metric_names,
  one_sample_summary(yearly_relative_effect),
  by = .(day_definition, metric)
]
scale_results[, `:=`(
  t_holm_p = p.adjust(t_p_value, method = "holm"),
  hac_holm_p = p.adjust(hac_p_value, method = "holm"),
  sign_holm_p = p.adjust(sign_p_value, method = "holm")
), by = day_definition]
scale_results[, scale_consistency := t_holm_p <= 0.025 &
                hac_holm_p <= 0.025 & sign_holm_p <= 0.025]

spatial_variance <- year_effects[
  analysis_role == "confirmatory" & metric == "spatial_variance",
  one_sample_summary(yearly_relative_effect),
  by = day_definition
]

descriptive_mean_summary <- function(values) {
  n <- length(values)
  estimate <- mean(values)
  standard_error <- sd(values) / sqrt(n)
  critical <- qt(0.975, df = n - 1L)
  data.table(
    years = n,
    estimate = estimate,
    standard_error = standard_error,
    ci_lower = estimate - critical * standard_error,
    ci_upper = estimate + critical * standard_error
  )
}

# Prespecified sensitivity definitions and thresholds. These reuse exactly the
# year-level estimator and are not part of the primary consistency decision.
sensitivity_specs <- list(
  threshold_70 = list(definition = "utc", lower = 0.25, upper = 0.70),
  threshold_80 = list(definition = "utc", lower = 0.25, upper = 0.80),
  stull_wbt = list(definition = "utc_stull", lower = 0.25, upper = 0.75),
  sitewise_daily_max = list(definition = "sitewise_max", lower = 0.25,
                            upper = 0.75)
)
sensitivity_years <- rbindlist(lapply(names(sensitivity_specs), function(label) {
  specification <- sensitivity_specs[[label]]
  if (!specification$definition %chin% unique(fields$day_definition)) return(NULL)
  fit <- one_definition(specification$definition, specification$lower,
                        specification$upper)
  copy(fit$profile)[, `:=`(
    sensitivity = label,
    lower_quantile = specification$lower,
    upper_quantile = specification$upper
  )]
}), use.names = TRUE, fill = TRUE)
sensitivity_results <- sensitivity_years[analysis_role == "confirmatory",
  one_sample_summary(yearly_profile_effect),
  by = .(sensitivity, day_definition, lower_quantile, upper_quantile)
]

# Continuous within-record association, summarized with the year as replicate.
continuous_years <- rbindlist(lapply(fits, function(fit) {
  slopes <- fit$metrics[, .(
    monthly_slope = cov(wbt_mean, profile_dispersion) / var(wbt_mean)
  ), by = .(file, year, analysis_role)]
  slopes[, .(
    yearly_slope = mean(monthly_slope),
    months = .N
  ), by = .(year, analysis_role)][, day_definition :=
    unique(fit$record$day_definition)]
}))
if (continuous_years[, any(months != 3L)]) stop("Incomplete continuous slopes")
continuous_results <- continuous_years[analysis_role == "confirmatory",
  descriptive_mean_summary(yearly_slope), by = day_definition]

# Early/late estimates are descriptive heterogeneity summaries, not a trend test.
heterogeneity_years <- copy(profile_years[analysis_role == "confirmatory"])
heterogeneity_years[, period := fifelse(year <= 2007L, "1991-2007",
                                         "2008-2025")]
heterogeneity_results <- heterogeneity_years[,
  descriptive_mean_summary(yearly_profile_effect),
  by = .(day_definition, period)]

# Descriptive associations with the retained land-surface variables. Each
# slope is estimated within a month, then averaged equally to the year level.
mechanism_years <- rbindlist(lapply(fits, function(fit) {
  if (is.null(fit$mechanisms)) return(NULL)
  joined <- merge(fit$metrics, fit$mechanisms, by = c("file", "date"))
  variables <- intersect(c("wind_speed_mean", "swvl1", "ssrd"), names(joined))
  rbindlist(lapply(variables, function(variable) {
    monthly <- joined[, .(
      monthly_slope = cov(profile_dispersion, get(variable)) / var(get(variable)),
      monthly_correlation = cor(profile_dispersion, get(variable))
    ), by = .(file, year, analysis_role)]
    monthly[, .(
      yearly_slope = mean(monthly_slope),
      yearly_correlation = mean(monthly_correlation),
      months = .N
    ), by = .(year, analysis_role)][, variable := variable]
  }))[, day_definition := unique(fit$record$day_definition)]
}), use.names = TRUE, fill = TRUE)
mechanism_results <- mechanism_years[analysis_role == "confirmatory", {
  slope_summary <- descriptive_mean_summary(yearly_slope)
  correlation_summary <- descriptive_mean_summary(yearly_correlation)
  data.table(
    years = slope_summary$years,
    mean_yearly_slope = slope_summary$estimate,
    slope_ci_lower = slope_summary$ci_lower,
    slope_ci_upper = slope_summary$ci_upper,
    mean_yearly_correlation = correlation_summary$estimate,
    correlation_ci_lower = correlation_summary$ci_lower,
    correlation_ci_upper = correlation_summary$ci_upper
  )
}, by = .(day_definition, variable)]

fwrite(record_effects, file.path(output_dir, "confirmatory_record_effects.csv"))
fwrite(year_effects, file.path(output_dir, "confirmatory_year_scale_effects.csv"))
fwrite(profile_years, file.path(output_dir, "confirmatory_year_profile_effects.csv"))
fwrite(primary, file.path(output_dir, "confirmatory_primary_results.csv"))
fwrite(scale_results, file.path(output_dir, "confirmatory_scale_results.csv"))
fwrite(spatial_variance,
       file.path(output_dir, "confirmatory_spatial_variance_results.csv"))
fwrite(rbindlist(lapply(fits, `[[`, "metadata"), idcol = "definition_index"),
       file.path(output_dir, "confirmatory_graph_metadata.csv"))
fwrite(sensitivity_years,
       file.path(output_dir, "confirmatory_sensitivity_year_effects.csv"))
fwrite(sensitivity_results,
       file.path(output_dir, "confirmatory_sensitivity_results.csv"))
fwrite(continuous_years,
       file.path(output_dir, "confirmatory_continuous_year_slopes.csv"))
fwrite(continuous_results,
       file.path(output_dir, "confirmatory_continuous_results.csv"))
fwrite(heterogeneity_results,
       file.path(output_dir, "confirmatory_period_results.csv"))
fwrite(mechanism_years,
       file.path(output_dir, "confirmatory_mechanism_year_effects.csv"))
fwrite(mechanism_results,
       file.path(output_dir, "confirmatory_mechanism_results.csv"))

cat("Confirmatory analysis complete.\n")
print(primary, digits = 4)
