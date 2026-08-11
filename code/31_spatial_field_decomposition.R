###############################################################################
# Spatial-field composites and an exact nodewise decomposition of the graph
# profile effect. The sum of the mapped node contributions equals the scalar
# profile estimate for each analysis role.
###############################################################################

library(data.table)
library(ggplot2)
library(maps)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))

input_dir <- file.path(project_dir, "data", "era5_confirmatory", "daily_fields")
output_dir <- file.path(project_dir, "output_confirmatory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = "^era5_land_[0-9]{4}_jja_daily_fields\\.csv\\.gz$",
  full.names = TRUE
)
if (length(files) != 35L) stop("Expected 35 yearly daily-field files")

keep <- c("year", "month", "record_id", "analysis_date", "day_definition",
          "site_id", "requested_lon", "requested_lat", "regional_mean_wbt",
          "wbt")
fields <- rbindlist(lapply(files, function(path) fread(path, select = keep)))
fields <- fields[day_definition == "utc", .(
  year, month, file = as.character(record_id),
  date = as.IDate(analysis_date), site_id,
  lon = requested_lon, lat = requested_lat,
  regional_mean_wbt, wbt
)]
setorder(fields, file, date, site_id)

if (fields[, uniqueN(site_id)] != 121L) stop("Expected 121 fixed sites")
counts <- fields[, .N, by = .(file, date)]
if (counts[, any(N != 121L)]) stop("Incomplete UTC spatial field")

daily_mean <- fields[, .(
  wbt_mean = unique(regional_mean_wbt),
  calculated_mean = mean(wbt)
), by = .(file, date)]
if (daily_mean[, max(abs(wbt_mean - calculated_mean))] > 1e-5) {
  stop("Stored and reconstructed regional means disagree")
}
classified <- classify_mean_quantiles(
  daily_mean[, .(file, date, wbt_mean)], lower = 0.25, upper = 0.75
)

sites <- unique(fields[, .(site_id, lon, lat)])[order(site_id)]
h_factors <- c(0.125, 0.25, 0.5, 1, 2)
operators <- make_graph_operators(sites, h_factors)
metric_names <- sprintf("graph_h_%s", format(h_factors, trim = TRUE))
scale_columns <- paste0(
  "contribution_h_",
  gsub("\\.", "_", format(h_factors, trim = TRUE))
)
discovery_years <- c(2015L, 2022L)

# Work record by record. Matrix operations calculate all daily node energies
# without materializing 121 x 121 pairwise arrays for every field.
records <- unique(fields[, .(file, year, month)])[order(year, month)]
record_maps <- rbindlist(lapply(seq_len(nrow(records)), function(k) {
  key <- records[k]
  x <- fields[file == key$file]
  dates <- sort(unique(x$date))
  y <- dcast(x, date ~ site_id, value.var = "wbt")
  if (!identical(y$date, dates)) stop("Date ordering failed in ", key$file)
  Y <- as.matrix(y[, -"date"])
  if (ncol(Y) != nrow(sites)) stop("Site ordering failed in ", key$file)

  regimes <- classified[file == key$file][match(dates, date), regime]
  if (anyNA(regimes)) stop("Missing regime in ", key$file)
  is_extreme <- regimes == "extreme"
  is_moderate <- regimes == "moderate"
  if (sum(is_extreme) < 2L || sum(is_moderate) < 2L) {
    stop("Insufficient regime size in ", key$file)
  }

  anomaly <- Y - rowMeans(Y)
  extreme_anomaly <- colMeans(anomaly[is_extreme, , drop = FALSE])
  moderate_anomaly <- colMeans(anomaly[is_moderate, , drop = FALSE])

  contribution_by_h <- vapply(operators, function(op) {
    degree <- rowSums(op$W)
    Y2 <- Y^2
    incident_energy <- sweep(Y2, 2, degree, "*") -
      2 * Y * (Y %*% op$W) + Y2 %*% op$W
    node_energy <- incident_energy / (4 * op$weight_sum)
    node_extreme <- colMeans(node_energy[is_extreme, , drop = FALSE])
    node_moderate <- colMeans(node_energy[is_moderate, , drop = FALSE])
    q_moderate <- sum(node_moderate)
    (node_extreme - node_moderate) / q_moderate
  }, numeric(nrow(sites)))
  colnames(contribution_by_h) <- scale_columns

  cbind(data.table(
    file = key$file,
    year = key$year,
    month = key$month,
    analysis_role = if (key$year %in% discovery_years) {
      "discovery"
    } else {
      "confirmatory"
    },
    sites,
    moderate_anomaly = moderate_anomaly,
    extreme_anomaly = extreme_anomaly,
    anomaly_difference = extreme_anomaly - moderate_anomaly,
    profile_contribution = rowMeans(contribution_by_h)
  ), as.data.table(contribution_by_h))
}))

# Preserve the primary estimator's weighting: equal months within a year,
# equal years within an analysis role, and equal graph bandwidths above.
map_columns <- c("moderate_anomaly", "extreme_anomaly",
                 "anomaly_difference", "profile_contribution",
                 scale_columns)
year_maps <- record_maps[, lapply(.SD, mean),
  by = .(analysis_role, year, site_id, lon, lat), .SDcols = map_columns]
role_maps <- year_maps[, lapply(.SD, mean),
  by = .(analysis_role, site_id, lon, lat), .SDcols = map_columns]

expected_file <- file.path(output_dir, "confirmatory_year_profile_effects.csv")
if (!file.exists(expected_file)) stop("Run 23_confirmatory_analysis.R first")
expected_years <- fread(expected_file)[day_definition == "utc"]
expected <- expected_years[, .(expected_effect = mean(yearly_profile_effect)),
                           by = analysis_role]
checks <- merge(
  role_maps[, .(mapped_effect = sum(profile_contribution)), by = analysis_role],
  expected, by = "analysis_role"
)
checks[, absolute_error := abs(mapped_effect - expected_effect)]
if (checks[, any(absolute_error > 1e-10)]) {
  stop("Node contributions do not recover the scalar profile estimates")
}

fwrite(role_maps,
       file.path(output_dir, "spatial_decomposition_role_maps.csv"))
fwrite(checks,
       file.path(output_dir, "spatial_decomposition_identity_check.csv"))

# The three smooth scale-attribution maps are display interpolations of these
# exact 121-node allocations.  Verify the raw-node identities before any
# smoothing so that the plotted node sums remain tied to the primary estimator.
scale_identity <- melt(
  role_maps[analysis_role == "confirmatory"],
  id.vars = c("analysis_role", "site_id", "lon", "lat"),
  measure.vars = scale_columns,
  variable.name = "scale_column", value.name = "node_contribution"
)[, .(mapped_effect = sum(node_contribution)),
  by = .(analysis_role, scale_column)]
scale_identity[, metric := metric_names[match(scale_column, scale_columns)]]
scale_reference <- fread(file.path(
  output_dir, "confirmatory_scale_results.csv"
))[day_definition == "utc" & metric %in% metric_names,
   .(metric, expected_effect = estimate)]
scale_identity <- merge(scale_identity, scale_reference, by = "metric")
scale_identity[, absolute_error := abs(mapped_effect - expected_effect)]
setcolorder(
  scale_identity,
  c("analysis_role", "metric", "scale_column", "mapped_effect",
    "expected_effect", "absolute_error")
)
if (nrow(scale_identity) != length(scale_columns) ||
    scale_identity[, any(absolute_error > 1e-10)]) {
  stop("Scale-specific node contributions do not recover primary estimates")
}
fwrite(
  scale_identity,
  file.path(output_dir, "spatial_decomposition_scale_identity_check.csv")
)

role_maps[, role_label := factor(
  analysis_role,
  levels = c("confirmatory", "discovery"),
  labels = c("Held-out confirmatory (33 summers)",
             "Discovery, excluded (2 summers)")
)]
world <- map_data("world")
world <- world[world$long >= 102 & world$long <= 128 &
                 world$lat >= 17 & world$lat <= 45, ]

theme_map <- theme_bw(base_size = 9.4) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 9.7, hjust = 0.5),
    axis.title = element_text(size = 8.8),
    axis.text = element_text(size = 7.7),
    legend.position = "bottom",
    legend.title = element_text(size = 8.2),
    legend.text = element_text(size = 7.5),
    legend.key.width = grid::unit(1.05, "cm"),
    plot.margin = margin(4, 2, 2, 2)
  )

map_panel <- function(data, variable, title, limit, legend_title,
                      show_x = TRUE, show_y = TRUE, sum_label = NULL) {
  p <- ggplot() +
    geom_polygon(
      data = world, aes(long, lat, group = group),
      fill = "#F2F2F2", colour = NA
    ) +
    geom_point(
      data = data, aes(lon, lat, fill = .data[[variable]]),
      shape = 22, size = 3.45, stroke = 0.18, colour = "#555555"
    ) +
    geom_path(
      data = world, aes(long, lat, group = group),
      colour = "#666666", linewidth = 0.25
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-limit, limit), name = legend_title,
      guide = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(4.0, "cm"),
        barheight = grid::unit(0.23, "cm")
      )
    ) +
    coord_quickmap(xlim = c(104.2, 125.0), ylim = c(19.3, 42.8),
                   expand = FALSE) +
    scale_x_continuous(breaks = c(105, 115, 125)) +
    scale_y_continuous(breaks = c(20, 30, 40)) +
    labs(
      x = if (show_x) "Longitude (degrees E)" else NULL,
      y = if (show_y) "Latitude (degrees N)" else NULL,
      title = title
    ) +
    theme_map
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  if (!is.null(sum_label)) {
    p <- p + annotate(
      "label", x = 104.7, y = 42.15, hjust = 0, vjust = 1,
      label = sum_label, size = 2.65, label.size = 0.18,
      fill = scales::alpha("white", 0.88)
    )
  }
  p
}

anomaly_limit <- max(abs(role_maps[, c(moderate_anomaly, extreme_anomaly)]))
difference_limit <- max(abs(role_maps$anomaly_difference))
contribution_limit <- 100 * max(abs(role_maps$profile_contribution))

titles <- c(
  "Moderate-day anomaly",
  "Extreme-day anomaly",
  "Extreme - moderate",
  "Profile contribution"
)

make_panels <- function(role) {
  d <- role_maps[analysis_role == role]
  effect <- 100 * sum(d$profile_contribution)
  d[, profile_contribution_pp := 100 * profile_contribution]
  list(
    map_panel(d, "moderate_anomaly", titles[1], anomaly_limit,
              "WBT anomaly (degrees C)", FALSE, TRUE),
    map_panel(d, "extreme_anomaly", titles[2], anomaly_limit,
              "WBT anomaly (degrees C)", FALSE, FALSE),
    map_panel(d, "anomaly_difference", titles[3], difference_limit,
              "Anomaly difference (degrees C)", TRUE, TRUE),
    map_panel(d, "profile_contribution_pp", titles[4], contribution_limit,
              "Contribution (percentage points)", TRUE, FALSE,
              sprintf("Sum = %+.2f%%", effect))
  )
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required for the spatial figure")
}
build_figure <- function(role, title, subtitle) {
  panels <- make_panels(role)
  assembled <- ((panels[[1]] | panels[[2]]) /
                  (panels[[3]] | panels[[4]]) +
                  patchwork::plot_layout(guides = "keep") +
                  patchwork::plot_annotation(
                    title = title,
                    subtitle = subtitle,
                    tag_levels = "a",
                    tag_prefix = "(",
                    tag_suffix = ")",
                    theme = theme(
                      plot.title = element_text(face = "bold", size = 12),
                      plot.subtitle = element_text(size = 9.2),
                      plot.tag = element_text(face = "bold", size = 10)
                    )
                  )) & theme(legend.position = "bottom")
  assembled
}

cat("Spatial identity checks:\n")
print(checks)
cat("Scale-specific spatial identity checks:\n")
print(scale_identity)
cat("Wrote spatial decomposition tables and identity checks.\n")
