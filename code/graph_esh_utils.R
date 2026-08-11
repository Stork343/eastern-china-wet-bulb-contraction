###############################################################################
# Multi-scale graph-dispersion tools for spatial homogenization analysis.
###############################################################################

library(data.table)

make_graph_operators <- function(sites,
                                 h_factors = c(0.125, 0.25, 0.5, 1, 2)) {
  sites <- copy(sites)[order(site_id)]
  coords <- project_coordinates_km(sites$lon, sites$lat)
  distance <- as.matrix(dist(coords))
  median_distance <- median(distance[lower.tri(distance)])

  lapply(h_factors, function(factor) {
    h <- factor * median_distance
    W <- exp(-(distance^2) / (2 * h^2))
    diag(W) <- 0
    L <- diag(rowSums(W)) - W
    weight_sum <- sum(W[upper.tri(W)])
    eigenvalues <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
    list(
      h_factor = factor,
      bandwidth_km = h,
      W = W,
      L = L,
      weight_sum = weight_sum,
      lambda_2 = sort(eigenvalues)[2],
      lambda_max = max(eigenvalues)
    )
  })
}

graph_dispersion <- function(y, operator) {
  as.numeric(crossprod(y, operator$L %*% y) /
               (2 * operator$weight_sum))
}

compute_graph_profile <- function(daily,
                                  h_factors = c(0.125, 0.25, 0.5, 1, 2)) {
  sites <- unique(daily[, .(site_id, lon, lat)])[order(site_id)]
  operators <- make_graph_operators(sites, h_factors)
  labels <- sprintf("graph_h_%s", format(h_factors, trim = TRUE))

  out <- daily[, {
    y <- wbt[order(site_id)]
    values <- vapply(operators, function(op) graph_dispersion(y, op), numeric(1))
    as.list(setNames(values, labels))
  }, by = .(file, date)]

  metadata <- rbindlist(lapply(seq_along(operators), function(i) {
    op <- operators[[i]]
    data.table(
      metric = labels[i],
      h_factor = op$h_factor,
      bandwidth_km = op$bandwidth_km,
      weight_sum = op$weight_sum,
      lambda_2 = op$lambda_2,
      lambda_max = op$lambda_max
    )
  }))
  list(metrics = out, metadata = metadata, operators = operators)
}

profile_effects <- function(classified, metric_names) {
  rbindlist(lapply(metric_names, function(metric) {
    ext <- classified[regime == "extreme", get(metric)]
    mod <- classified[regime == "moderate", get(metric)]
    data.table(
      metric = metric,
      mean_extreme = mean(ext),
      mean_moderate = mean(mod),
      delta = mean(ext) - mean(mod),
      relative_delta = mean(ext) / mean(mod) - 1
    )
  }))
}

shift_profile_null <- function(classified, metric_names, B, seed) {
  strata <- lapply(unique(classified$file), function(f) {
    x <- classified[file == f][order(date)]
    list(
      regime = x$regime,
      outcomes = as.matrix(x[, ..metric_names]),
      n = nrow(x)
    )
  })

  set.seed(seed)
  null <- matrix(NA_real_, nrow = B, ncol = length(metric_names),
                 dimnames = list(NULL, metric_names))
  for (b in seq_len(B)) {
    ext_sum <- mod_sum <- numeric(length(metric_names))
    ext_n <- mod_n <- 0L
    for (x in strata) {
      offset <- sample.int(x$n, 1L) - 1L
      idx <- ((seq_len(x$n) - 1L + offset) %% x$n) + 1L
      shifted <- x$outcomes[idx, , drop = FALSE]
      is_ext <- x$regime == "extreme"
      is_mod <- x$regime == "moderate"
      ext_sum <- ext_sum + colSums(shifted[is_ext, , drop = FALSE])
      mod_sum <- mod_sum + colSums(shifted[is_mod, , drop = FALSE])
      ext_n <- ext_n + sum(is_ext)
      mod_n <- mod_n + sum(is_mod)
    }
    ext_mean <- ext_sum / ext_n
    mod_mean <- mod_sum / mod_n
    null[b, ] <- ext_mean / mod_mean - 1
  }
  null
}

circular_shift_profile_test <- function(metrics, metric_names,
                                        B = 9999L, seed = 20260870L,
                                        lower = 0.25, upper = 0.75) {
  classified <- classify_mean_quantiles(metrics, lower, upper)
  observed_table <- profile_effects(classified, metric_names)
  observed <- observed_table$relative_delta
  null <- shift_profile_null(classified, metric_names, B, seed)

  null_center <- colMeans(null)
  null_scale <- apply(null, 2, sd)
  observed_z <- (observed - null_center) / null_scale
  null_z <- sweep(sweep(null, 2, null_center), 2, null_scale, "/")

  # Directed omnibus alternative: the mean relative dispersion change across
  # prespecified bandwidths is negative. This statistic is fixed before the
  # randomization distribution is generated, as required by the group test.
  observed_omnibus <- mean(observed)
  null_omnibus <- rowMeans(null)
  omnibus_p <- (1 + sum(null_omnibus <= observed_omnibus)) / (B + 1)

  # Single-step minimum-relative-effect adjustment respects the common
  # percentage scale and dependence across bandwidths.
  null_min <- apply(null, 1, min)
  adjusted_p <- vapply(observed, function(effect) {
    (1 + sum(null_min <= effect)) / (B + 1)
  }, numeric(1))

  observed_table[, `:=`(
    z_score = observed_z,
    unadjusted_p = vapply(seq_along(observed), function(k) {
      (1 + sum(null[, k] <= observed[k])) / (B + 1)
    }, numeric(1)),
    fwer_p = adjusted_p
  )]

  list(
    effects = observed_table,
    omnibus_statistic = observed_omnibus,
    omnibus_p = omnibus_p,
    null_relative_effects = null,
    null_omnibus = null_omnibus,
    B = B,
    classified = classified
  )
}
