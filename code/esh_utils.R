###############################################################################
# Shared utilities for the corrected ESH analysis.
###############################################################################

library(data.table)

load_consistent_era5 <- function(data_dir) {
  files <- list.files(data_dir, pattern = "^era5_[0-9]{6}_consistent\\.csv$",
                      full.names = TRUE)
  if (length(files) != 6L) stop("Expected six consistent ERA5 files in ", data_dir)

  out <- rbindlist(lapply(files, function(f) {
    x <- fread(f)
    x[, file := basename(f)]
    x
  }), use.names = TRUE, fill = TRUE)

  out[, time_utc := as.POSIXct(time, tz = "UTC")]
  if (out[, anyNA(time_utc)]) stop("Failed to parse UTC timestamps")
  setorder(out, file, time_utc, site_id)

  site_check <- out[, .(n_coord = uniqueN(paste(lon, lat))), by = site_id]
  if (site_check[, any(n_coord != 1L)]) stop("site_id is not a stable spatial key")
  if (out[, uniqueN(paste(lon, lat))] != 121L) stop("Expected 121 fixed sites")
  out
}

aggregate_daily_composite_max <- function(hourly, wbt_column = "wbt") {
  hourly[, date := as.IDate(time_utc, tz = "UTC")]
  daily <- hourly[, .(
    wbt = max(get(wbt_column), na.rm = TRUE),
    t2m = max(t2m, na.rm = TRUE)
  ), by = .(file, date, site_id, lon, lat)]

  counts <- daily[, .N, by = .(file, date)]
  if (counts[, any(N != 121L)]) stop("Incomplete daily spatial field detected")
  setorder(daily, file, date, site_id)
  daily
}

aggregate_daily_peak_fields <- function(hourly, wbt_column = "wbt") {
  hourly[, date := as.IDate(time_utc, tz = "UTC")]
  regional <- hourly[, .(regional_mean_wbt = mean(get(wbt_column))),
                     by = .(file, date, time_utc)]
  setorder(regional, file, date, -regional_mean_wbt, time_utc)
  peaks <- regional[, .SD[1L], by = .(file, date)]
  daily <- merge(
    hourly,
    peaks[, .(file, date, time_utc, regional_mean_wbt)],
    by = c("file", "date", "time_utc")
  )[, .(file, date, peak_time_utc = time_utc, site_id, lon, lat,
        wbt = get(wbt_column), t2m, regional_mean_wbt)]

  counts <- daily[, .N, by = .(file, date)]
  if (counts[, any(N != 121L)]) stop("Incomplete synchronous daily field")
  setorder(daily, file, date, site_id)
  daily
}

# Backward-compatible alias; corrected analyses use synchronous peak fields.
aggregate_daily_fields <- aggregate_daily_peak_fields

project_coordinates_km <- function(lon, lat) {
  lat0 <- mean(lat) * pi / 180
  cbind(x = lon * 111.32 * cos(lat0), y = lat * 110.57)
}

compute_field_metrics <- function(daily, h_factor = 0.5) {
  sites <- unique(daily[, .(site_id, lon, lat)])[order(site_id)]
  if (nrow(sites) != 121L) stop("Expected 121 sites")
  coords <- project_coordinates_km(sites$lon, sites$lat)
  distance <- as.matrix(dist(coords))
  bandwidth_km <- median(distance[lower.tri(distance)]) * h_factor
  weights <- exp(-(distance^2) / (2 * bandwidth_km^2))

  metrics <- daily[, {
    y <- wbt[order(site_id)]
    if (length(y) != nrow(sites)) stop("Incomplete field")
    centered <- y - mean(y)
    weighted_matrix <- tcrossprod(centered) * weights / length(y)
    weighted_matrix <- (weighted_matrix + t(weighted_matrix)) / 2
    lambda_weighted <- eigen(weighted_matrix, symmetric = TRUE,
                             only.values = TRUE)$values[1]

    .(
      n_sites = length(y),
      wbt_mean = mean(y),
      spatial_variance = mean(centered^2),
      lambda_weighted = lambda_weighted
    )
  }, by = .(file, date)]

  attr(metrics, "bandwidth_km") <- bandwidth_km
  setorder(metrics, file, date)
  metrics
}

classify_mean_quantiles <- function(metrics, lower = 0.25, upper = 0.75) {
  out <- copy(metrics)
  out[, c("q_lower", "q_upper") := {
    q <- quantile(wbt_mean, c(lower, upper), names = FALSE, type = 7)
    .(q[1], q[2])
  }, by = file]
  out[, regime := fcase(
    wbt_mean >= q_upper, "extreme",
    wbt_mean >= q_lower & wbt_mean < q_upper, "moderate",
    default = "excluded_low"
  )]
  out
}

estimate_delta <- function(classified, metric = "spatial_variance") {
  ext <- classified[regime == "extreme", get(metric)]
  mod <- classified[regime == "moderate", get(metric)]
  if (length(ext) < 2L || length(mod) < 2L) return(NA_real_)
  mean(ext) - mean(mod)
}

circular_block_indices <- function(n, block_length) {
  n_blocks <- ceiling(n / block_length)
  starts <- sample.int(n, n_blocks, replace = TRUE)
  idx <- unlist(lapply(starts, function(start) {
    ((start - 1L + seq_len(block_length) - 1L) %% n) + 1L
  }), use.names = FALSE)
  idx[seq_len(n)]
}

stratified_block_resample <- function(metrics, block_length) {
  pieces <- lapply(unique(metrics$file), function(f) {
    x <- metrics[file == f][order(date)]
    x[circular_block_indices(nrow(x), block_length)]
  })
  rbindlist(pieces, use.names = TRUE)
}

block_bootstrap_delta <- function(metrics, metric = "spatial_variance",
                                  block_length = 7L, B = 1999L,
                                  seed = 20260802L,
                                  lower = 0.25, upper = 0.75) {
  set.seed(seed)
  observed_class <- classify_mean_quantiles(metrics, lower, upper)
  observed <- estimate_delta(observed_class, metric)
  if (!is.finite(observed)) stop("Observed delta is not finite")

  boot <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    sample_b <- stratified_block_resample(metrics, block_length)
    boot[b] <- estimate_delta(
      classify_mean_quantiles(sample_b, lower, upper), metric
    )
  }
  boot <- boot[is.finite(boot)]
  if (length(boot) < 0.95 * B) stop("Too many invalid bootstrap replicates")

  q <- quantile(boot, c(0.025, 0.975), names = FALSE)
  null_statistics <- boot - observed
  p_lower <- (1 + sum(null_statistics <= observed)) / (length(boot) + 1)

  list(
    observed = observed,
    basic_ci = c(lower = 2 * observed - q[2], upper = 2 * observed - q[1]),
    percentile_ci = c(lower = q[1], upper = q[2]),
    p_value = p_lower,
    bootstrap_se = sd(boot),
    boot = boot,
    classified = observed_class,
    block_length = block_length,
    B = length(boot)
  )
}

monthly_effects <- function(metrics, metric = "spatial_variance") {
  classified <- classify_mean_quantiles(metrics)
  classified[, .(
    n_extreme = sum(regime == "extreme"),
    n_moderate = sum(regime == "moderate"),
    mean_extreme = mean(get(metric)[regime == "extreme"]),
    mean_moderate = mean(get(metric)[regime == "moderate"]),
    delta = mean(get(metric)[regime == "extreme"]) -
      mean(get(metric)[regime == "moderate"])
  ), by = file]
}

circular_shift_test <- function(metrics, metric = "spatial_variance",
                                B = 9999L, seed = 20260850L,
                                lower = 0.25, upper = 0.75) {
  classified <- classify_mean_quantiles(metrics, lower, upper)
  observed <- estimate_delta(classified, metric)
  files <- unique(classified$file)
  strata <- lapply(files, function(f) {
    x <- classified[file == f][order(date)]
    list(regime = x$regime, outcome = x[[metric]], n = nrow(x))
  })

  set.seed(seed)
  null <- numeric(B)
  for (b in seq_len(B)) {
    ext_sum <- mod_sum <- 0
    ext_n <- mod_n <- 0L
    for (x in strata) {
      offset <- sample.int(x$n, 1L) - 1L
      idx <- ((seq_len(x$n) - 1L + offset) %% x$n) + 1L
      shifted_outcome <- x$outcome[idx]
      is_ext <- x$regime == "extreme"
      is_mod <- x$regime == "moderate"
      ext_sum <- ext_sum + sum(shifted_outcome[is_ext])
      mod_sum <- mod_sum + sum(shifted_outcome[is_mod])
      ext_n <- ext_n + sum(is_ext)
      mod_n <- mod_n + sum(is_mod)
    }
    null[b] <- ext_sum / ext_n - mod_sum / mod_n
  }

  list(
    observed = observed,
    p_value = (1 + sum(null <= observed)) / (B + 1),
    null = null,
    B = B,
    classified = classified
  )
}

record_level_inference <- function(metrics, metric = "spatial_variance",
                                   conf_level = 0.95,
                                   lower = 0.25, upper = 0.75) {
  classified <- classify_mean_quantiles(metrics, lower, upper)
  effects <- classified[, .(
    delta = mean(get(metric)[regime == "extreme"]) -
      mean(get(metric)[regime == "moderate"])
  ), by = file]$delta
  n <- length(effects)
  estimate <- mean(effects)
  se <- sd(effects) / sqrt(n)
  alpha <- 1 - conf_level
  critical <- qt(1 - alpha / 2, df = n - 1L)
  list(
    estimate = estimate,
    ci = c(lower = estimate - critical * se,
           upper = estimate + critical * se),
    p_value_lower = pt(estimate / se, df = n - 1L),
    n_records = n,
    effects = effects
  )
}

circular_shift_slope_test <- function(metrics, metric = "spatial_variance",
                                      B = 9999L, seed = 20260860L) {
  x <- copy(metrics)
  x[, mean_centered := wbt_mean - mean(wbt_mean), by = file]
  observed <- coef(lm(x[[metric]] ~ x$mean_centered + factor(x$file)))[2]
  strata <- split(x[order(file, date)], by = "file", keep.by = TRUE)

  set.seed(seed)
  null <- numeric(B)
  for (b in seq_len(B)) {
    shifted <- rbindlist(lapply(strata, function(z) {
      n <- nrow(z)
      offset <- sample.int(n, 1L) - 1L
      idx <- ((seq_len(n) - 1L + offset) %% n) + 1L
      out <- copy(z)
      out[, (metric) := z[[metric]][idx]]
      out
    }))
    null[b] <- coef(lm(shifted[[metric]] ~ shifted$mean_centered +
                         factor(shifted$file)))[2]
  }

  list(
    slope = unname(observed),
    p_value = (1 + sum(null <= observed)) / (B + 1),
    null = null,
    B = B
  )
}
