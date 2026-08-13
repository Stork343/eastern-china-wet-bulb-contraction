###############################################################################
# Publication maps from the primary 121-node analysis.
#
# The thin-plate REML predictions in this script are display coordinates only.
# Every anomaly, node allocation, estimand, identity check, and annotation is
# read from or checked against the primary 121-node confirmatory outputs.
# Predictions are never substituted for nodes in estimation or inference.
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

required_packages <- c("patchwork", "sp", "jsonlite", "ragg")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

confirmatory_dir <- file.path(project_dir, "results")
output_dir <- file.path(project_dir, "results")
manuscript_figure_dir <- file.path(project_dir, "manuscript", "figures")
preview_dir <- file.path(manuscript_figure_dir, "png")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  spatial_nodes = file.path(
    confirmatory_dir, "spatial_decomposition_role_maps.csv"
  ),
  profile_identity = file.path(
    confirmatory_dir, "spatial_decomposition_identity_check.csv"
  ),
  scale_identity = file.path(
    confirmatory_dir, "spatial_decomposition_scale_identity_check.csv"
  ),
  scale_results = file.path(
    confirmatory_dir, "confirmatory_scale_results.csv"
  ),
  annual_scale_effects = file.path(
    confirmatory_dir, "confirmatory_year_scale_effects.csv"
  ),
  graph_metadata = file.path(
    confirmatory_dir, "confirmatory_graph_metadata.csv"
  ),
  primary_manifest = file.path(
    project_dir, "data", "grid", "eastern_china_121_sites.csv"
  )
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Missing primary-map input files: ", paste(missing_inputs, collapse = ", "))
}

all_nodes <- fread(input_files[["spatial_nodes"]])
nodes <- all_nodes[analysis_role == "confirmatory"]
setorder(nodes, site_id)
required_node_columns <- c(
  "site_id", "lon", "lat", "moderate_anomaly", "extreme_anomaly",
  "anomaly_difference", "profile_contribution"
)
if (!all(required_node_columns %in% names(nodes))) {
  stop("Primary spatial-decomposition output lacks required columns")
}
if (nrow(nodes) != 121L || uniqueN(nodes$site_id) != 121L ||
    anyNA(nodes[, ..required_node_columns])) {
  stop("Expected 121 complete, unique primary analysis nodes")
}

# The manifest validates site identity and coordinates; it supplies no plotted
# outcome or allocation value.
manifest <- fread(input_files[["primary_manifest"]])
setorder(manifest, site_id)
if (nrow(manifest) != 121L ||
    !identical(nodes$site_id, manifest$site_id) ||
    max(abs(nodes$lon - manifest$lon)) > 1e-12 ||
    max(abs(nodes$lat - manifest$lat)) > 1e-12) {
  stop("Spatial output does not match the frozen 121-site manifest")
}

metadata <- fread(input_files[["graph_metadata"]])[
  definition_index == 1L,
  .(metric, h_factor, bandwidth_km)
]
metadata <- unique(metadata)
setorder(metadata, h_factor)
if (nrow(metadata) != 5L) stop("Expected five primary graph scales")
scale_columns <- paste0(
  "contribution_h_",
  gsub("\\.", "_", format(metadata$h_factor, trim = TRUE))
)
if (!all(scale_columns %in% names(nodes))) {
  stop("Run 31_spatial_field_decomposition.R to retain scale allocations")
}

profile_identity <- fread(input_files[["profile_identity"]])[
  analysis_role == "confirmatory"
]
if (nrow(profile_identity) != 1L ||
    profile_identity$absolute_error > 1e-10) {
  stop("Primary profile node identity is unavailable or failed")
}
profile_node_sum <- sum(nodes$profile_contribution)
if (abs(profile_node_sum - profile_identity$expected_effect) > 1e-10) {
  stop("Primary profile allocation no longer recovers the scalar estimand")
}

scale_identity <- fread(input_files[["scale_identity"]])[
  analysis_role == "confirmatory"
]
scale_reference <- fread(input_files[["scale_results"]])[
  day_definition == "utc" & metric %in% metadata$metric,
  .(metric, estimate)
]
scale_check <- merge(
  metadata[, .(metric, h_factor, bandwidth_km, scale_column = scale_columns)],
  scale_identity[, .(metric, mapped_effect, expected_effect, absolute_error)],
  by = "metric"
)
scale_check <- merge(scale_check, scale_reference, by = "metric")
setorder(scale_check, h_factor)
if (nrow(scale_check) != 5L ||
    scale_check[, any(absolute_error > 1e-10)] ||
    scale_check[, any(abs(mapped_effect - estimate) > 1e-10)] ||
    scale_check[, any(abs(expected_effect - estimate) > 1e-10)]) {
  stop("Scale-specific primary node identities failed")
}

###############################################################################
# Supported high-resolution display grid.
###############################################################################

display_resolution_degrees <- 0.1
nearest_node_gate_km <- 140
spline_k <- 30L
lat0 <- mean(nodes$lat) * pi / 180

project_xy <- function(lon, lat) {
  data.table(
    x = lon * 111.32 * cos(lat0),
    y = lat * 110.57
  )
}

prediction_candidates <- CJ(
  lon = seq(min(nodes$lon), max(nodes$lon), by = display_resolution_degrees),
  lat = seq(min(nodes$lat), max(nodes$lat), by = display_resolution_degrees)
)
candidate_xy <- project_xy(
  prediction_candidates$lon, prediction_candidates$lat
)
node_xy <- project_xy(nodes$lon, nodes$lat)
nearest_km <- rep(Inf, nrow(prediction_candidates))
for (k in seq_len(nrow(node_xy))) {
  distance <- sqrt(
    (candidate_xy$x - node_xy$x[k])^2 +
      (candidate_xy$y - node_xy$y[k])^2
  )
  nearest_km <- pmin(nearest_km, distance)
}

hull_indices <- chull(nodes$lon, nodes$lat)
hull <- nodes[c(hull_indices, hull_indices[1L]), .(lon, lat)]
inside_hull <- sp::point.in.polygon(
  prediction_candidates$lon, prediction_candidates$lat,
  hull$lon, hull$lat
) > 0L
on_land <- !is.na(map.where(
  "world", prediction_candidates$lon, prediction_candidates$lat
))
prediction_grid <- prediction_candidates[
  inside_hull & on_land & nearest_km <= nearest_node_gate_km
]
prediction_xy <- project_xy(prediction_grid$lon, prediction_grid$lat)
prediction_grid[, `:=`(x = prediction_xy$x, y = prediction_xy$y)]
if (nrow(prediction_grid) < 10000L) {
  stop("Supported display grid is unexpectedly sparse")
}

surface_audit <- list()
fit_surface <- function(node_data, variable, figure_id, panel_id,
                        reference_sum = NA_real_) {
  xy <- project_xy(node_data$lon, node_data$lat)
  model_data <- data.table(
    value = node_data[[variable]],
    x = xy$x,
    y = xy$y
  )
  fit <- gam(
    value ~ s(x, y, bs = "tp", k = spline_k),
    data = model_data,
    method = "REML",
    select = TRUE
  )
  surface <- copy(prediction_grid)
  surface[, value := as.numeric(predict(fit, newdata = surface))]
  smooth_edf <- unname(summary(fit)$s.table[1L, "edf"])
  raw_sum <- sum(model_data$value)
  audit_row <- data.table(
    figure = figure_id,
    panel = panel_id,
    source_variable = variable,
    analysis_nodes = nrow(model_data),
    prediction_cells = nrow(surface),
    display_resolution_degrees = display_resolution_degrees,
    spline_basis = "thin-plate regression spline",
    spline_k = spline_k,
    smoothing_method = "REML",
    term_selection = TRUE,
    smooth_edf = smooth_edf,
    raw_node_min = min(model_data$value),
    raw_node_max = max(model_data$value),
    predicted_min = min(surface$value),
    predicted_max = max(surface$value),
    raw_node_sum = raw_sum,
    reference_sum = reference_sum,
    raw_identity_error = if (is.finite(reference_sum)) {
      abs(raw_sum - reference_sum)
    } else {
      NA_real_
    },
    land_mask = TRUE,
    convex_hull_clip = TRUE,
    nearest_node_gate_km = nearest_node_gate_km,
    display_only = TRUE,
    predictions_enter_estimation_or_inference = FALSE
  )
  surface_audit[[length(surface_audit) + 1L]] <<- audit_row
  list(surface = surface, edf = smooth_edf)
}

###############################################################################
# Shared publication styling and geographic base.
###############################################################################

ink <- "#24303A"
blue <- "#1F5A85"
orange <- "#D55E00"
paper <- "#FFFFFF"
grid_colour <- "#D9DEE3"
land_colour <- "#EDF0F2"

world <- map_data("world")
world <- world[
  world$long >= 102 & world$long <= 128 &
    world$lat >= 17 & world$lat <= 45,
]

map_theme <- theme_bw(base_size = 9.4, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    panel.background = element_rect(fill = paper, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = grid_colour, linewidth = 0.20),
    panel.border = element_rect(colour = "#8D98A0", linewidth = 0.32),
    plot.title = element_text(
      face = "bold", size = 9.8, hjust = 0.5, margin = margin(b = 3)
    ),
    plot.tag = element_text(face = "bold", size = 10.2),
    axis.title = element_text(size = 8.8),
    axis.text = element_text(size = 7.7, colour = ink),
    legend.position = "bottom",
    legend.title = element_text(size = 8.0),
    legend.text = element_text(size = 7.4),
    legend.margin = margin(t = -2),
    plot.margin = margin(4, 3, 2, 3)
  )

base_surface_map <- function(surface, node_data, title, tag,
                             legend_title, colour_limit,
                             show_x = TRUE, show_y = TRUE,
                             sum_label = NULL) {
  panel <- ggplot() +
    geom_polygon(
      data = world, aes(long, lat, group = group),
      fill = land_colour, colour = NA
    ) +
    geom_raster(
      data = surface, aes(lon, lat, fill = value), interpolate = TRUE
    ) +
    geom_path(
      data = world, aes(long, lat, group = group),
      colour = "#5E6870", linewidth = 0.28
    ) +
    geom_point(
      data = node_data, aes(lon, lat), inherit.aes = FALSE,
      shape = 16, size = 0.46,
      colour = "#30383E", alpha = 0.64
    ) +
    scale_fill_gradient2(
      low = blue, mid = paper, high = orange, midpoint = 0,
      limits = c(-colour_limit, colour_limit), oob = scales::squish,
      name = legend_title,
      guide = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(4.0, "cm"),
        barheight = grid::unit(0.23, "cm")
      )
    ) +
    coord_quickmap(
      xlim = c(104.2, 125.1), ylim = c(19.3, 42.8), expand = FALSE
    ) +
    scale_x_continuous(breaks = c(105, 110, 115, 120, 125)) +
    scale_y_continuous(breaks = c(20, 25, 30, 35, 40)) +
    labs(
      x = if (show_x) "Longitude (degrees E)" else NULL,
      y = if (show_y) "Latitude (degrees N)" else NULL,
      title = title,
      tag = tag
    ) +
    map_theme
  if (!show_x) {
    panel <- panel + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank()
    )
  }
  if (!show_y) {
    panel <- panel + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
  }
  if (!is.null(sum_label)) {
    panel <- panel + annotate(
      "label", x = 104.65, y = 42.15, hjust = 0, vjust = 1,
      label = sum_label, size = 2.65, label.size = 0.16,
      label.padding = grid::unit(0.12, "lines"),
      fill = scales::alpha(paper, 0.90), colour = ink
    )
  }
  panel
}

###############################################################################
# Four-panel spatial decomposition.
###############################################################################

nodes[, profile_contribution_pp := 100 * profile_contribution]
composite_specs <- list(
  list(
    variable = "moderate_anomaly", panel = "Middle-day spatially centred WBT",
    legend = "WBT minus regional mean (degrees C)", reference = NA_real_
  ),
  list(
    variable = "extreme_anomaly", panel = "High-day spatially centred WBT",
    legend = "WBT minus regional mean (degrees C)", reference = NA_real_
  ),
  list(
    variable = "anomaly_difference", panel = "High - middle",
    legend = "Anomaly difference (degrees C)", reference = NA_real_
  ),
  list(
    variable = "profile_contribution_pp", panel = "Profile allocation",
    legend = "Node allocation (percentage points)",
    reference = 100 * profile_identity$expected_effect
  )
)
composite_surfaces <- lapply(composite_specs, function(spec) {
  fit_surface(
    nodes, spec$variable, "primary_spatial_decomposition", spec$panel,
    spec$reference
  )
})
anomaly_limit <- max(abs(c(
  nodes$moderate_anomaly, nodes$extreme_anomaly,
  composite_surfaces[[1L]]$surface$value,
  composite_surfaces[[2L]]$surface$value
)))
difference_limit <- max(abs(c(
  nodes$anomaly_difference,
  composite_surfaces[[3L]]$surface$value
)))
allocation_limit <- max(abs(c(
  nodes$profile_contribution_pp,
  composite_surfaces[[4L]]$surface$value
)))
composite_limits <- c(
  anomaly_limit, anomaly_limit, difference_limit, allocation_limit
)
composite_panels <- lapply(seq_along(composite_specs), function(k) {
  spec <- composite_specs[[k]]
  base_surface_map(
    composite_surfaces[[k]]$surface, nodes,
    spec$panel, sprintf("(%s)", LETTERS[k]), spec$legend,
    composite_limits[k],
    show_x = k >= 3L,
    show_y = k %in% c(1L, 3L),
    sum_label = if (k == 4L) {
      sprintf("Raw-node sum = %+.2f%%", 100 * profile_node_sum)
    } else {
      NULL
    }
  )
})
primary_spatial_figure <- (
  (composite_panels[[1L]] | composite_panels[[2L]]) /
    (composite_panels[[3L]] | composite_panels[[4L]])
) & theme(legend.position = "bottom")

###############################################################################
# Exact scale-specific allocations and annual evidence.
###############################################################################

scale_indices <- c(1L, 3L, 5L)
scale_map_data <- lapply(scale_indices, function(index) {
  column <- scale_columns[index]
  d <- copy(nodes)
  d[, node_allocation_pp := 100 * get(column)]
  surface <- fit_surface(
    d, "node_allocation_pp", "primary_scale_attribution",
    sprintf("%s-km graph", format(round(metadata$bandwidth_km[index]), big.mark = ",")),
    100 * scale_check$estimate[index]
  )$surface
  list(nodes = d, surface = surface, index = index)
})
scale_limit <- max(abs(unlist(lapply(scale_map_data, function(item) {
  c(item$nodes$node_allocation_pp, item$surface$value)
}))))
scale_panels <- lapply(seq_along(scale_map_data), function(k) {
  item <- scale_map_data[[k]]
  index <- item$index
  base_surface_map(
    item$surface, item$nodes,
    sprintf("%s-km graph", format(round(metadata$bandwidth_km[index]), big.mark = ",")),
    sprintf("(%s)", LETTERS[k]),
    "Node allocation (percentage points)", scale_limit,
    show_x = TRUE, show_y = k == 1L,
    sum_label = sprintf(
      "Raw-node sum = %+.2f%%",
      sum(item$nodes$node_allocation_pp)
    )
  )
})

annual_scale <- fread(input_files[["annual_scale_effects"]])[
  day_definition == "utc" & analysis_role == "confirmatory" &
    metric %in% metadata$metric
]
annual_scale <- merge(
  annual_scale,
  metadata[, .(metric, bandwidth_km)],
  by = "metric"
)
annual_scale[, bandwidth_label := factor(
  sprintf("%s km", format(round(bandwidth_km), big.mark = ",")),
  levels = sprintf(
    "%s km",
    format(round(sort(unique(metadata$bandwidth_km))), big.mark = ",")
  )
)]
annual_summary <- fread(input_files[["scale_results"]])[
  day_definition == "utc" & metric %in% metadata$metric
]
annual_summary <- merge(
  annual_summary,
  metadata[, .(metric, bandwidth_km)],
  by = "metric"
)
annual_summary[, bandwidth_label := factor(
  sprintf("%s km", format(round(bandwidth_km), big.mark = ",")),
  levels = levels(annual_scale$bandwidth_label)
)]
if (nrow(annual_scale) != 33L * 5L ||
    uniqueN(annual_scale$year) != 33L || nrow(annual_summary) != 5L) {
  stop("Held-out annual scale evidence is incomplete")
}

distribution_theme <- theme_bw(base_size = 9.4, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.26),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(colour = "#8D98A0", linewidth = 0.32),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.tag = element_text(face = "bold", size = 10.2),
    axis.title = element_text(size = 8.9),
    axis.text = element_text(size = 7.8, colour = ink),
    legend.position = "none",
    plot.margin = margin(4, 4, 3, 4)
  )
if (exists("figure3", inherits = FALSE) && inherits(figure3, "ggplot")) {
  # When called by 38_jrssc_main_figures.R, reuse its already-audited annual
  # panel verbatim apart from the redundant standalone title and subtitle.
  distribution_panel_source <- "in-memory figure3 from 38_jrssc_main_figures.R"
  distribution_panel <- figure3 +
    labs(title = NULL, subtitle = NULL, tag = "(D)") +
    theme(
      plot.tag = element_text(face = "bold", size = 10.2),
      plot.margin = margin(4, 4, 3, 4)
    )
} else {
  # Standalone execution reproduces the same panel from the frozen 121-node
  # annual outputs so that this script remains independently runnable.
  distribution_panel_source <- "reconstructed from confirmatory annual output"
  distribution_panel <- ggplot() +
    annotate(
      "rect", xmin = -40, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = "#F3F7FA", colour = NA
    ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#68727A") +
    geom_point(
      data = annual_scale,
      aes(100 * yearly_relative_effect, bandwidth_label),
      position = position_jitter(width = 0, height = 0.13, seed = 20260809),
      shape = 21, size = 1.8, stroke = 0.25, colour = paper,
      fill = "#8FB6D1", alpha = 0.72
    ) +
    geom_errorbarh(
      data = annual_summary,
      aes(y = bandwidth_label, xmin = 100 * ci_lower, xmax = 100 * ci_upper),
      height = 0, linewidth = 1.0, colour = blue
    ) +
    geom_point(
      data = annual_summary,
      aes(100 * estimate, bandwidth_label),
      shape = 21, size = 3.5, stroke = 0.80, colour = blue, fill = paper
    ) +
    geom_text(
      data = annual_summary,
      aes(
        x = 11.5, y = bandwidth_label,
        label = sprintf("%d/33", negative_years)
      ),
      hjust = 0.5, size = 3.0, colour = ink
    ) +
    scale_x_continuous(
      limits = c(-40, 16), breaks = seq(-40, 10, 10),
      labels = function(x) sprintf("%g%%", x)
    ) +
    labs(
      x = "Annual high-to-middle graph-dispersion contrast",
      y = "Gaussian graph bandwidth",
      title = NULL,
      subtitle = NULL,
      tag = "(D)"
    ) +
    distribution_theme
}

primary_scale_figure <- (
  (scale_panels[[1L]] | scale_panels[[2L]] | scale_panels[[3L]]) /
    distribution_panel +
    patchwork::plot_layout(heights = c(1.02, 0.68), guides = "collect")
) & theme(legend.position = "bottom")

###############################################################################
# Stable vector outputs, manuscript copies, previews, and audit metadata.
###############################################################################

save_pdf <- function(path, figure, width, height) {
  ggsave(
    path, figure, width = width, height = height,
    device = cairo_pdf, bg = paper
  )
}
save_png <- function(path, figure, width, height) {
  ggsave(
    path, figure, width = width, height = height, units = "in",
    device = ragg::agg_png, dpi = 320, bg = paper
  )
}

project_relative <- function(path) {
  normalised <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- paste0(
    normalizePath(project_dir, winslash = "/", mustWork = TRUE), "/"
  )
  ifelse(
    startsWith(normalised, root),
    substring(normalised, nchar(root) + 1L),
    normalised
  )
}

output_files <- c(
  primary_spatial_pdf = file.path(
    output_dir, "fig_primary_spatial_decomposition.pdf"
  ),
  primary_scale_pdf = file.path(
    output_dir, "fig_primary_scale_attribution.pdf"
  ),
  main_multiscale_pdf = file.path(
    output_dir, "fig3_multiscale_evidence.pdf"
  ),
  manuscript_multiscale_pdf = file.path(
    manuscript_figure_dir, "figure03_multiscale_evidence.pdf"
  ),
  manuscript_spatial_pdf = file.path(
    manuscript_figure_dir, "figure04_primary_spatial_decomposition.pdf"
  ),
  spatial_preview = file.path(
    preview_dir, "figure04_primary_spatial_decomposition.png"
  ),
  multiscale_preview = file.path(
    preview_dir, "figure03_multiscale_evidence.png"
  )
)
save_pdf(output_files[["primary_spatial_pdf"]], primary_spatial_figure, 8.3, 7.8)
save_pdf(output_files[["main_multiscale_pdf"]], primary_scale_figure, 8.3, 7.4)
portable_copies <- c(
  primary_scale_pdf = file.copy(
    output_files[["main_multiscale_pdf"]],
    output_files[["primary_scale_pdf"]], overwrite = TRUE
  ),
  manuscript_multiscale_pdf = file.copy(
    output_files[["main_multiscale_pdf"]],
    output_files[["manuscript_multiscale_pdf"]], overwrite = TRUE
  ),
  manuscript_spatial_pdf = file.copy(
    output_files[["primary_spatial_pdf"]],
    output_files[["manuscript_spatial_pdf"]], overwrite = TRUE
  )
)
if (!all(portable_copies)) stop("Failed to copy canonical map PDFs")
save_png(output_files[["spatial_preview"]], primary_spatial_figure, 8.3, 7.8)
save_png(output_files[["multiscale_preview"]], primary_scale_figure, 8.3, 7.4)

surface_audit_table <- rbindlist(surface_audit, fill = TRUE)
if (nrow(surface_audit_table) != 7L ||
    surface_audit_table[, any(analysis_nodes != 121L)] ||
    surface_audit_table[
      is.finite(reference_sum), any(raw_identity_error > 1e-10)
    ] ||
    surface_audit_table[, any(predictions_enter_estimation_or_inference)]) {
  stop("Primary smooth-map audit failed")
}
surface_audit_file <- file.path(
  output_dir, "primary_smooth_map_surface_audit.csv"
)
fwrite(surface_audit_table, surface_audit_file)

build_audit <- list(
  analysis_support = "primary 121-node analysis grid",
  analysis_role = "primary analysis (33 evaluation summers)",
  display_interpolation = list(
    method = "low-rank thin-plate regression spline",
    basis = "tp",
    k = spline_k,
    smoothing_parameter = "REML",
    term_selection = TRUE,
    grid_resolution_degrees = display_resolution_degrees,
    support_mask = paste(
      "maps::world land, within the 121-node convex hull, and no farther than",
      paste0(nearest_node_gate_km, " km from a node")
    ),
    display_only = TRUE,
    predictions_enter_estimation_or_inference = FALSE
  ),
  exact_raw_node_identities = list(
    profile_percent = 100 * profile_node_sum,
    profile_absolute_error = abs(
      profile_node_sum - profile_identity$expected_effect
    ),
    selected_scale_percent = setNames(
      as.list(100 * scale_check$mapped_effect[scale_indices]),
      sprintf("%.0f_km", metadata$bandwidth_km[scale_indices])
    ),
    selected_scale_max_absolute_error = max(
      scale_check$absolute_error[scale_indices]
    )
  ),
  annual_distribution = list(
    summers = uniqueN(annual_scale$year),
    graph_scales = uniqueN(annual_scale$metric),
    panel_source = distribution_panel_source,
    interpolation_used = FALSE
  ),
  numeric_value_sources = as.list(project_relative(
    input_files[names(input_files) != "primary_manifest"]
  )),
  coordinate_validation_source = project_relative(
    input_files[["primary_manifest"]]
  ),
  prohibited_dense_inputs_used = FALSE,
  output_files = as.list(project_relative(output_files)),
  surface_audit_csv = project_relative(surface_audit_file)
)
build_audit_file <- file.path(output_dir, "primary_smooth_map_build_audit.json")
jsonlite::write_json(build_audit, build_audit_file, pretty = TRUE, auto_unbox = TRUE)

cat("Primary raw-node profile sum (%):\n")
print(100 * profile_node_sum)
cat("Primary selected-scale raw-node sums (%):\n")
print(scale_check[h_factor %in% metadata$h_factor[scale_indices], .(
  bandwidth_km, mapped_percent = 100 * mapped_effect,
  identity_error = absolute_error
)])
cat("Surface prediction and EDF audit:\n")
print(surface_audit_table[, .(
  figure, panel, smooth_edf, raw_node_min, raw_node_max,
  predicted_min, predicted_max, raw_node_sum, raw_identity_error
)])
cat("Wrote primary-only smooth spatial figures and audit metadata.\n")
