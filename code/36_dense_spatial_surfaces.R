###############################################################################
# Exact dense-grid node decomposition and descriptive smooth spatial surfaces.
#
# Inference remains on the 465 sampled nodes. Thin-plate REML surfaces are
# clipped to land within 140 km of a sampled site and are used only to display
# the empirical spatial pattern at publication scale.
###############################################################################

library(data.table)
library(ggplot2)
library(maps)
library(mgcv)

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
          "label_mean_wbt", "wbt")
fields <- rbindlist(lapply(files, function(path) fread(path, select = keep)))
fields[, `:=`(
  file = as.character(record_id), date = as.IDate(analysis_date),
  lon = requested_lon, lat = requested_lat
)]
setorder(fields, analysis_definition, file, date, site_id)
if (fields[, uniqueN(site_id)] != 465L) stop("Expected 465 dense sites")

manifest <- fread(manifest_file)
sites <- manifest[, .(site_id = dense_site_id, lon, lat)][order(site_id)]
metadata <- fread(file.path(
  primary_output, "confirmatory_graph_metadata.csv"
))[definition_index == 1L]
metadata <- unique(metadata[, .(metric, h_factor, bandwidth_km)])
setorder(metadata, h_factor)

make_fixed_operators <- function(sites, metadata) {
  coords <- project_coordinates_km(sites$lon, sites$lat)
  distance <- as.matrix(dist(coords))
  lapply(seq_len(nrow(metadata)), function(k) {
    h <- metadata$bandwidth_km[k]
    W <- exp(-(distance^2) / (2 * h^2))
    diag(W) <- 0
    list(metric = metadata$metric[k], bandwidth_km = h, W = W,
         weight_sum = sum(W[upper.tri(W)]))
  })
}
operators <- make_fixed_operators(sites, metadata)
scale_columns <- paste0(
  "contribution_h_",
  gsub("\\.", "_", format(metadata$h_factor, trim = TRUE))
)

definition_config <- data.table(
  analysis_definition = c("primary_grid_peak", "dense_grid_peak"),
  configuration = c("dense_465_fixed_labels", "dense_465_recomputed")
)

one_definition <- function(definition, configuration) {
  x_all <- fields[analysis_definition == definition]
  records <- unique(x_all[, .(file, year, month, analysis_role)])[
    order(year, month)
  ]
  rbindlist(lapply(seq_len(nrow(records)), function(k) {
    key <- records[k]
    x <- x_all[file == key$file]
    dates <- sort(unique(x$date))
    labels <- unique(x[, .(date, label_mean_wbt)])[order(date)]
    if (!identical(labels$date, dates)) stop("Label date mismatch")
    q <- quantile(labels$label_mean_wbt, c(0.25, 0.75),
                  names = FALSE, type = 7)
    regimes <- fcase(
      labels$label_mean_wbt >= q[2], "extreme",
      labels$label_mean_wbt >= q[1] & labels$label_mean_wbt < q[2],
      "moderate",
      default = "excluded_low"
    )
    is_extreme <- regimes == "extreme"
    is_moderate <- regimes == "moderate"

    wide <- dcast(x[, .(date, site_id, wbt)], date ~ site_id,
                  value.var = "wbt")
    value_columns <- as.character(sites$site_id)
    Y <- as.matrix(wide[, ..value_columns])
    if (!identical(wide$date, dates) || ncol(Y) != 465L) {
      stop("Dense record matrix mismatch")
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
      (node_extreme - node_moderate) / sum(node_moderate)
    }, numeric(nrow(sites)))
    colnames(contribution_by_h) <- scale_columns

    cbind(data.table(
      configuration = configuration,
      file = key$file, year = key$year, month = key$month,
      analysis_role = key$analysis_role,
      sites,
      moderate_anomaly = moderate_anomaly,
      extreme_anomaly = extreme_anomaly,
      anomaly_difference = extreme_anomaly - moderate_anomaly,
      profile_contribution = rowMeans(contribution_by_h)
    ), as.data.table(contribution_by_h))
  }))
}

record_maps <- rbindlist(lapply(seq_len(nrow(definition_config)), function(k) {
  one_definition(definition_config$analysis_definition[k],
                 definition_config$configuration[k])
}))
map_columns <- c("moderate_anomaly", "extreme_anomaly",
                 "anomaly_difference", "profile_contribution",
                 scale_columns)
year_maps <- record_maps[, lapply(.SD, mean),
  by = .(configuration, analysis_role, year, site_id, lon, lat),
  .SDcols = map_columns]
role_maps <- year_maps[, lapply(.SD, mean),
  by = .(configuration, analysis_role, site_id, lon, lat),
  .SDcols = map_columns]

expected <- fread(file.path(output_dir, "dense_year_profile_effects.csv"))[,
  .(expected_effect = mean(yearly_profile_effect)),
  by = .(configuration, analysis_role)]
identity <- merge(
  role_maps[, .(mapped_effect = sum(profile_contribution)),
            by = .(configuration, analysis_role)],
  expected, by = c("configuration", "analysis_role")
)
identity[, absolute_error := abs(mapped_effect - expected_effect)]
if (identity[, any(absolute_error > 1e-10)]) {
  stop("Dense node contributions do not recover the scalar estimates")
}
fwrite(role_maps, file.path(output_dir, "dense_spatial_role_maps.csv"))
fwrite(identity, file.path(output_dir, "dense_spatial_identity_check.csv"))

scale_identity <- melt(
  role_maps[configuration == "dense_465_fixed_labels" &
              analysis_role == "confirmatory"],
  id.vars = c("configuration", "analysis_role", "site_id", "lon", "lat"),
  measure.vars = scale_columns,
  variable.name = "scale_column", value.name = "node_contribution"
)[, .(mapped_effect = sum(node_contribution)),
  by = .(configuration, analysis_role, scale_column)]
scale_identity[, metric := metadata$metric[match(scale_column, scale_columns)]]
scale_reference <- fread(file.path(output_dir, "dense_scale_results.csv"))[
  configuration == "dense_465_fixed_labels", .(metric, expected_effect = estimate)
]
scale_identity <- merge(scale_identity, scale_reference, by = "metric")
scale_identity[, absolute_error := abs(mapped_effect - expected_effect)]
if (scale_identity[, any(absolute_error > 1e-10)]) {
  stop("Scale-specific node contributions do not recover dense estimates")
}
fwrite(scale_identity,
       file.path(output_dir, "dense_scale_identity_check.csv"))

# Prediction support: a 0.1-degree display grid, land only, and no farther than
# 140 km from a sampled site. The distance gate prevents unsupported coastal
# or corner extrapolation.
lat0 <- mean(sites$lat) * pi / 180
project_xy <- function(lon, lat) {
  data.table(x = lon * 111.32 * cos(lat0), y = lat * 110.57)
}
prediction_grid <- CJ(
  lon = seq(105.0, 125.0, by = 0.1),
  lat = seq(20.4, 42.0, by = 0.1)
)
grid_xy <- project_xy(prediction_grid$lon, prediction_grid$lat)
site_xy <- project_xy(sites$lon, sites$lat)
nearest_km <- rep(Inf, nrow(prediction_grid))
for (k in seq_len(nrow(site_xy))) {
  distance <- sqrt((grid_xy$x - site_xy$x[k])^2 +
                     (grid_xy$y - site_xy$y[k])^2)
  nearest_km <- pmin(nearest_km, distance)
}
prediction_grid[, `:=`(
  x = grid_xy$x,
  y = grid_xy$y,
  nearest_site_km = nearest_km,
  land = !is.na(map.where("world", lon, lat))
)]
prediction_grid <- prediction_grid[land == TRUE & nearest_site_km <= 140]

smooth_one <- function(site_data, variable) {
  xy <- project_xy(site_data$lon, site_data$lat)
  model_data <- data.table(
    value = site_data[[variable]], x = xy$x, y = xy$y
  )
  fit <- gam(value ~ s(x, y, bs = "tp", k = 80),
             data = model_data, method = "REML", select = TRUE)
  surface <- copy(prediction_grid)
  surface[, value := as.numeric(predict(fit, newdata = surface))]
  attr(surface, "edf") <- sum(fit$edf)
  surface
}

world <- map_data("world")
world <- world[world$long >= 102 & world$long <= 128 &
                 world$lat >= 17 & world$lat <= 45, ]
theme_map <- theme_bw(base_size = 9.6) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.8),
    legend.position = "bottom",
    legend.title = element_text(size = 8.4),
    legend.text = element_text(size = 7.6),
    plot.margin = margin(4, 3, 2, 3)
  )

surface_panel <- function(site_data, variable, title, legend_title,
                          show_x, show_y, sum_label = NULL) {
  surface <- smooth_one(site_data, variable)
  limit <- max(abs(c(surface$value, site_data[[variable]])))
  p <- ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "#EEEEEE", colour = NA) +
    geom_raster(data = surface, aes(lon, lat, fill = value), interpolate = TRUE) +
    geom_path(data = world, aes(long, lat, group = group),
              colour = "#555555", linewidth = 0.25) +
    geom_point(data = site_data, aes(lon, lat), inherit.aes = FALSE,
               colour = "#222222", size = 0.28, alpha = 0.55) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-limit, limit), name = legend_title,
      guide = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(4.0, "cm"),
        barheight = grid::unit(0.23, "cm")
      )
    ) +
    coord_quickmap(xlim = c(104.2, 125.1), ylim = c(19.3, 42.8),
                   expand = FALSE) +
    scale_x_continuous(breaks = c(105, 115, 125)) +
    scale_y_continuous(breaks = c(20, 30, 40)) +
    labs(
      x = if (show_x) "Longitude (degrees E)" else NULL,
      y = if (show_y) "Latitude (degrees N)" else NULL,
      title = title
    ) + theme_map
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  if (!is.null(sum_label)) {
    p <- p + annotate(
      "label", x = 104.7, y = 42.15, hjust = 0, vjust = 1,
      label = sum_label, size = 2.7, label.size = 0.18,
      fill = scales::alpha("white", 0.88)
    )
  }
  list(plot = p, surface = surface)
}

make_surface_figure <- function(configuration, role, filename, title_suffix) {
  config_value <- configuration
  role_value <- role
  d <- role_maps[
    configuration == config_value & analysis_role == role_value
  ]
  if (nrow(d) != 465L) stop("Missing dense role map")
  d[, profile_contribution_pp := 100 * profile_contribution]
  effect <- sum(d$profile_contribution_pp)
  panels <- list(
    surface_panel(d, "moderate_anomaly", "Middle-day anomaly",
                  "WBT anomaly (degrees C)", FALSE, TRUE),
    surface_panel(d, "extreme_anomaly", "High-day anomaly",
                  "WBT anomaly (degrees C)", FALSE, FALSE),
    surface_panel(d, "anomaly_difference", "High - middle",
                  "Anomaly difference (degrees C)", TRUE, TRUE),
    surface_panel(d, "profile_contribution_pp", "Profile contribution",
                  "Node contribution (percentage points)", TRUE, FALSE,
                  sprintf("Node sum = %+.2f%%", effect))
  )
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required")
  }
  figure <- (((panels[[1]]$plot | panels[[2]]$plot) /
              (panels[[3]]$plot | panels[[4]]$plot) +
              patchwork::plot_layout(guides = "keep") +
              patchwork::plot_annotation(
                tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
                caption = paste(
                  title_suffix,
                  "Thin-plate REML surfaces are descriptive; dots mark sampled sites."
                ),
                theme = theme(
                  plot.tag = element_text(face = "bold", size = 10),
                  plot.caption = element_text(size = 8, hjust = 0)
                )
              ))) & theme(legend.position = "bottom")
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), figure,
         width = 7.8, height = 8.0, device = cairo_pdf)
  edf <- vapply(panels, function(item) attr(item$surface, "edf"), numeric(1))
  data.table(configuration, analysis_role = role,
             variable = c("moderate_anomaly", "extreme_anomaly",
                          "anomaly_difference", "profile_contribution_pp"),
             smooth_edf = edf, node_effect = effect)
}

smooth_audit <- rbind(
  make_surface_figure(
    "dense_465_fixed_labels", "confirmatory",
    "fig10_dense_fixed_smooth_surfaces",
    "465-site support with primary peak times and labels."
  )
)
fwrite(smooth_audit, file.path(output_dir, "dense_surface_smoothing_audit.csv"))

###############################################################################
# Scale-specific spatial attribution on one common colour scale.
###############################################################################

scale_indices <- c(1L, 3L, 5L)
scale_data <- role_maps[
  configuration == "dense_465_fixed_labels" &
    analysis_role == "confirmatory"
]
if (nrow(scale_data) != 465L) stop("Missing scale-attribution map data")

scale_surfaces <- lapply(scale_columns[scale_indices], function(variable) {
  d <- copy(scale_data)
  d[, value_pp := 100 * get(variable)]
  list(nodes = d, surface = smooth_one(d, "value_pp"))
})
scale_limit <- max(abs(unlist(lapply(scale_surfaces, function(item) {
  c(item$nodes$value_pp, item$surface$value)
}))))

scale_map_panel <- function(item, index, show_x, show_y, tag) {
  h <- metadata$bandwidth_km[index]
  node_sum <- sum(item$nodes$value_pp)
  ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "#EEEEEE", colour = NA) +
    geom_raster(data = item$surface, aes(lon, lat, fill = value),
                interpolate = TRUE) +
    geom_path(data = world, aes(long, lat, group = group),
              colour = "#555555", linewidth = 0.25) +
    geom_point(data = item$nodes, aes(lon, lat), inherit.aes = FALSE,
               colour = "#222222", size = 0.27, alpha = 0.5) +
    annotate("label", x = 104.7, y = 42.15, hjust = 0, vjust = 1,
             label = sprintf("Node sum = %+.2f%%", node_sum), size = 2.65,
             label.size = 0.18, fill = scales::alpha("white", 0.88)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-scale_limit, scale_limit),
      name = "Node contribution (percentage points)",
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5,
                              barwidth = grid::unit(4.4, "cm"),
                              barheight = grid::unit(0.23, "cm"))
    ) +
    coord_quickmap(xlim = c(104.2, 125.1), ylim = c(19.3, 42.8),
                   expand = FALSE) +
    scale_x_continuous(breaks = c(105, 110, 115, 120, 125)) +
    scale_y_continuous(breaks = c(20, 25, 30, 35, 40)) +
    labs(x = if (show_x) "Longitude (degrees E)" else NULL,
         y = if (show_y) "Latitude (degrees N)" else NULL,
         title = sprintf("%.0f-km graph", h), tag = tag) +
    theme_map +
    theme(
      plot.tag = element_text(face = "bold", size = 10),
      plot.tag.position = c(0.01, 0.99),
      axis.text.x = if (show_x) element_text(size = 7.8) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.text.y = if (show_y) element_text(size = 7.8) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank()
    )
}

# Longitude-averaged anomaly profiles per latitude band, plotted as anomaly
# against latitude so that the steep topographic gradient near 33-34 degrees
# north appears as a monotone step rather than as a cusp.
latitude_profile <- scale_data[, .(
  moderate_anomaly = mean(moderate_anomaly),
  extreme_anomaly = mean(extreme_anomaly)
), by = lat]
latitude_profile <- melt(latitude_profile, id.vars = "lat",
                         variable.name = "regime", value.name = "anomaly")
latitude_profile[, regime := factor(
  regime, levels = c("moderate_anomaly", "extreme_anomaly"),
  labels = c("Middle days", "High days")
)]
setorder(latitude_profile, regime, lat)
latitude_profile[, smooth_anomaly := as.numeric(predict(
  loess(anomaly ~ lat, span = 0.55, degree = 1),
  newdata = data.frame(lat = lat)
)), by = regime]
p_latitude <- ggplot(latitude_profile, aes(lat, anomaly, colour = regime)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#6B747B") +
  geom_point(size = 1.1, alpha = 0.48) +
  geom_line(aes(y = smooth_anomaly), linewidth = 0.82, lineend = "round") +
  scale_colour_manual(values = c("Middle days" = "#D55E00",
                                 "High days" = "#1F5A85"), name = "Regime") +
  scale_x_continuous(breaks = c(20, 25, 30, 35, 40),
                     limits = c(20.2, 42.2)) +
  labs(x = "Latitude (degrees N)",
       y = "Longitude-averaged WBT anomaly (degrees C)",
       title = "Latitudinal anomaly profile", tag = "(D)") +
  theme_map + theme(plot.tag = element_text(face = "bold", size = 10),
                    plot.tag.position = c(0.045, 0.96))

scale_panels <- lapply(seq_along(scale_indices), function(k) {
  scale_map_panel(scale_surfaces[[k]], scale_indices[k],
                  show_x = TRUE, show_y = k == 1L,
                  tag = sprintf("(%s)", LETTERS[k]))
})
scale_figure <- ((scale_panels[[1]] | scale_panels[[2]] | scale_panels[[3]]) /
                   p_latitude) +
  patchwork::plot_layout(guides = "collect", heights = c(1.0, 0.68)) +
  patchwork::plot_annotation(
    caption = paste(
      "465-site held-out analysis with primary peak times and labels.",
      "All three maps use one symmetric colour scale; dots mark sampled sites."
    ),
    theme = theme(plot.caption = element_text(size = 8, hjust = 0))
  ) & theme(legend.position = "bottom")
ggsave(file.path(output_dir, "fig11_dense_scale_attribution.pdf"),
       scale_figure, width = 8.3, height = 7.4, device = cairo_pdf)

scale_smooth_audit <- data.table(
  metric = metadata$metric[scale_indices],
  bandwidth_km = metadata$bandwidth_km[scale_indices],
  smooth_edf = vapply(scale_surfaces,
                      function(item) attr(item$surface, "edf"), numeric(1)),
  mapped_effect = vapply(scale_surfaces,
                         function(item) sum(item$nodes$value_pp) / 100,
                         numeric(1)),
  common_colour_limit_pp = scale_limit
)
fwrite(scale_smooth_audit,
       file.path(output_dir, "dense_scale_smoothing_audit.csv"))

cat("Dense spatial identity checks:\n")
print(identity)
cat("Wrote dense composite and scale-attribution surface figures.\n")
