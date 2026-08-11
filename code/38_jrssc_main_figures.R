###############################################################################
# Unified main-text figures for the JRSS C manuscript.
#
# Each main-text figure now has one statistical message and one coordinate
# system.  The shared colour-blind-safe palette and uncertainty conventions
# remain, but secondary diagnostics are left to tables and text.
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

corrected_dir <- file.path(project_dir, "output_corrected")
confirmatory_dir <- file.path(project_dir, "output_confirmatory")
dense_dir <- file.path(project_dir, "output_dense")
noaa_extension_dir <- file.path(project_dir, "output_noaa_extension")
output_dir <- file.path(project_dir, "output_jrssc")
manuscript_figure_dir <- file.path(project_dir, "manuscript", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_figure_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#24303A"
blue <- "#1F5A85"
light_blue <- "#A9C5DA"
orange <- "#D55E00"
teal <- "#008C72"
gold <- "#E6A700"
purple <- "#76528A"
paper <- "#FFFFFF"
grid <- "#D9DEE3"
muted <- "#8B98A3"
pale <- "#DCE5EC"

jr_theme <- function(base_size = 9.4) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = base_size + 0.6,
                                margin = margin(b = 5)),
      plot.title.position = "plot",
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = ink),
      panel.grid.major.y = element_line(colour = grid, linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "#8D98A0", linewidth = 0.32),
      axis.ticks = element_line(colour = "#8D98A0", linewidth = 0.32),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", colour = ink),
      plot.margin = margin(5, 6, 5, 6)
    )
}

save_figure <- function(name, figure, width, height) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")), figure,
         width = width, height = height, device = cairo_pdf,
         bg = paper)
}

###############################################################################
# Main-text spatial figure: the primary 121-site held-out decomposition.
###############################################################################

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required for the spatial figure")
}

spatial_maps <- fread(file.path(
  confirmatory_dir, "spatial_decomposition_role_maps.csv"
))[analysis_role == "confirmatory"]
if (nrow(spatial_maps) != 121L || uniqueN(spatial_maps$site_id) != 121L) {
  stop("Expected 121 sites in the confirmatory spatial decomposition")
}
spatial_maps[, node_allocation_pp := 100 * profile_contribution]
node_sum <- spatial_maps[, sum(node_allocation_pp)]
if (round(node_sum, 2) != -7.28) {
  stop("The primary node allocation does not sum to -7.28%")
}

map_background <- map_data("world")
map_background <- map_background[
  map_background$long >= 102 & map_background$long <= 128 &
    map_background$lat >= 17 & map_background$lat <= 45,
]

map_theme <- theme_bw(base_size = 9.2, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = grid, linewidth = 0.20),
    plot.title = element_text(face = "bold", size = 9.7, hjust = 0.5,
                              margin = margin(b = 3)),
    axis.title = element_text(size = 8.7),
    axis.text = element_text(size = 7.7, colour = ink),
    legend.position = "bottom",
    legend.title = element_text(size = 7.9),
    legend.text = element_text(size = 7.3),
    legend.key.width = grid::unit(1.02, "cm"),
    legend.margin = margin(t = -2),
    plot.margin = margin(4, 3, 2, 3)
  )

primary_map_panel <- function(data, variable, title, limit, legend_title,
                              show_x = TRUE, show_y = TRUE,
                              sum_label = NULL) {
  panel <- ggplot() +
    geom_polygon(
      data = map_background, aes(long, lat, group = group),
      fill = "#F1F3F4", colour = NA
    ) +
    geom_point(
      data = data, aes(lon, lat, fill = .data[[variable]]),
      shape = 22, size = 3.45, stroke = 0.20, colour = "#4F5961"
    ) +
    geom_path(
      data = map_background, aes(long, lat, group = group),
      colour = "#68727A", linewidth = 0.28
    ) +
    scale_fill_gradient2(
      low = blue, mid = paper, high = orange, midpoint = 0,
      limits = c(-limit, limit), oob = scales::squish,
      name = legend_title,
      guide = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(3.7, "cm"),
        barheight = grid::unit(0.22, "cm")
      )
    ) +
    coord_quickmap(
      xlim = c(104.2, 125.0), ylim = c(19.3, 42.8), expand = FALSE
    ) +
    scale_x_continuous(breaks = c(105, 115, 125)) +
    scale_y_continuous(breaks = c(20, 30, 40)) +
    labs(
      x = if (show_x) "Longitude (degrees E)" else NULL,
      y = if (show_y) "Latitude (degrees N)" else NULL,
      title = title
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
      "label", x = 104.7, y = 42.15, hjust = 0, vjust = 1,
      label = sum_label, size = 2.75, label.size = 0.18,
      fill = scales::alpha(paper, 0.90), colour = ink
    )
  }
  panel
}

anomaly_limit <- spatial_maps[, max(abs(c(
  moderate_anomaly, extreme_anomaly
)))]
difference_limit <- spatial_maps[, max(abs(anomaly_difference))]
allocation_limit <- spatial_maps[, max(abs(node_allocation_pp))]

spatial_panels <- list(
  primary_map_panel(
    spatial_maps, "moderate_anomaly", "Middle-day anomaly",
    anomaly_limit, "WBT anomaly (degrees C)", show_x = FALSE
  ),
  primary_map_panel(
    spatial_maps, "extreme_anomaly", "High-day anomaly",
    anomaly_limit, "WBT anomaly (degrees C)",
    show_x = FALSE, show_y = FALSE
  ),
  primary_map_panel(
    spatial_maps, "anomaly_difference", "High minus middle",
    difference_limit, "Anomaly difference (degrees C)"
  ),
  primary_map_panel(
    spatial_maps, "node_allocation_pp", "Node allocation",
    allocation_limit, "Allocation (percentage points)", show_y = FALSE,
    sum_label = sprintf("Node sum = %.2f%%", node_sum)
  )
)

spatial_figure <- (
  (spatial_panels[[1]] | spatial_panels[[2]]) /
    (spatial_panels[[3]] | spatial_panels[[4]]) +
    patchwork::plot_layout(guides = "keep") +
    patchwork::plot_annotation(
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 10, colour = ink)
      )
    )
) & theme(legend.position = "bottom")
save_figure("fig_primary_spatial_decomposition", spatial_figure, 8.3, 7.4)

###############################################################################
# Main-text energy decomposition across the five prespecified graph scales.
###############################################################################

energy <- fread(file.path(
  confirmatory_dir, "extended_energy_decomposition.csv"
))
if (nrow(energy) != 5L || max(energy$identity_error) > 1e-10) {
  stop("Energy-decomposition output is incomplete or fails its identity")
}
energy_long <- rbindlist(list(
  energy[, .(
    bandwidth_km, component = "Total raw-field contrast",
    estimate = total_effect_estimate,
    ci_lower = total_effect_ci_lower,
    ci_upper = total_effect_ci_upper
  )],
  energy[, .(
    bandwidth_km, component = "Anomaly-energy component",
    estimate = anomaly_energy_component_estimate,
    ci_lower = anomaly_energy_component_ci_lower,
    ci_upper = anomaly_energy_component_ci_upper
  )],
  energy[, .(
    bandwidth_km, component = "Climatology-anomaly cross component",
    estimate = climatology_anomaly_cross_component_estimate,
    ci_lower = climatology_anomaly_cross_component_ci_lower,
    ci_upper = climatology_anomaly_cross_component_ci_upper
  )]
))
component_levels <- c(
  "Total raw-field contrast",
  "Anomaly-energy component",
  "Climatology-anomaly cross component"
)
energy_long[, component := factor(component, levels = component_levels)]
energy_long[, bandwidth_label := factor(
  sprintf("%s", format(round(bandwidth_km), big.mark = ",")),
  levels = sprintf(
    "%s", format(round(sort(unique(bandwidth_km))), big.mark = ",")
  )
)]

energy_figure <- ggplot(
  energy_long,
  aes(bandwidth_label, 100 * estimate, colour = component,
      shape = component, linetype = component, group = component)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#6B747B") +
  geom_errorbar(
    aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper),
    width = 0.10, linewidth = 0.72, alpha = 0.82
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 3.1, stroke = 0.75, fill = paper) +
  scale_colour_manual(
    values = c(
      "Total raw-field contrast" = ink,
      "Anomaly-energy component" = teal,
      "Climatology-anomaly cross component" = orange
    ), name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Total raw-field contrast" = 21,
      "Anomaly-energy component" = 24,
      "Climatology-anomaly cross component" = 22
    ), name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Total raw-field contrast" = "solid",
      "Anomaly-energy component" = "dashed",
      "Climatology-anomaly cross component" = "dotdash"
    ), name = NULL
  ) +
  scale_y_continuous(labels = function(x) sprintf("%g%%", x)) +
  labs(
    x = "Gaussian graph bandwidth (km)",
    y = "High-minus-middle contrast component"
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.25),
    legend.position = "none"
  )
save_figure("fig_energy_decomposition", energy_figure, 8.3, 4.9)

###############################################################################
# Supplementary bandwidth-resolution and spatial-convergence diagnostics.
###############################################################################

dense_profile <- fread(file.path(
  confirmatory_dir, "extended_dense_bandwidth_profile.csv"
))
dense_curve <- dense_profile[source == "31_point_curve"]
fixed_points <- dense_profile[source == "prespecified_scale"]
if (nrow(dense_curve) != 31L || nrow(fixed_points) != 5L) {
  stop("Expected 31 dense bandwidths and five prespecified scales")
}
bandwidth_breaks <- sort(fixed_points$bandwidth_km)

dense_bandwidth_figure <- ggplot(
  dense_curve, aes(bandwidth_km, 100 * estimate)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#6B747B") +
  geom_ribbon(
    aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper),
    fill = light_blue, alpha = 0.48, colour = NA
  ) +
  geom_line(colour = blue, linewidth = 1.05) +
  geom_point(
    data = fixed_points,
    aes(bandwidth_km, 100 * estimate), inherit.aes = FALSE,
    shape = 21, size = 3.6, stroke = 1.0, colour = orange, fill = paper
  ) +
  scale_x_log10(
    breaks = bandwidth_breaks,
    labels = function(x) format(round(x), big.mark = ",", scientific = FALSE)
  ) +
  scale_y_continuous(labels = function(x) sprintf("%g%%", x)) +
  labs(
    x = "Gaussian graph bandwidth (km; log scale)",
    y = "Mean high-minus-middle graph-dispersion contrast"
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.25),
    legend.position = "none"
  )
save_figure("supp_dense_bandwidth_profile", dense_bandwidth_figure, 8.3, 4.8)

convergence <- fread(file.path(
  dense_dir, "extended_spatial_convergence.csv"
))
expected_configurations <- c(
  "primary_121", "latitude_refined", "longitude_refined", "dense_465"
)
if (!setequal(unique(convergence$configuration), expected_configurations) ||
    convergence[, uniqueN(bandwidth_key)] != 5L) {
  stop("Spatial-convergence output has unexpected configurations or scales")
}
convergence_summary <- convergence[, .(
  sites = unique(sites),
  max_difference_pp = 100 * max(abs(difference_from_dense))
), by = configuration]
convergence_summary[, configuration_label := fcase(
  configuration == "primary_121",
    sprintf("Primary (%d sites; max |difference| %.2f pp)",
            sites, max_difference_pp),
  configuration == "latitude_refined",
    sprintf("Latitude-refined (%d; %.2f pp)", sites, max_difference_pp),
  configuration == "longitude_refined",
    sprintf("Longitude-refined (%d; %.2f pp)", sites, max_difference_pp),
  configuration == "dense_465",
    sprintf("Dense reference (%d sites)", sites)
)]
convergence <- merge(
  convergence,
  convergence_summary[, .(configuration, configuration_label)],
  by = "configuration", all.x = TRUE, sort = FALSE
)
configuration_levels <- convergence_summary[
  match(expected_configurations, configuration), configuration_label
]
convergence[, configuration_label := factor(
  configuration_label, levels = configuration_levels
)]
convergence[, bandwidth_label := factor(
  format(round(bandwidth_key), big.mark = ",", scientific = FALSE),
  levels = format(
    round(sort(unique(bandwidth_key))), big.mark = ",", scientific = FALSE
  )
)]
convergence_colours <- setNames(
  c(orange, purple, teal, ink), configuration_levels
)
convergence_shapes <- setNames(c(21, 24, 22, 23), configuration_levels)
convergence_linetypes <- setNames(
  c("solid", "dashed", "dotdash", "longdash"), configuration_levels
)

convergence_figure <- ggplot(
  convergence,
  aes(bandwidth_label, 100 * difference_from_dense,
      colour = configuration_label, shape = configuration_label,
      linetype = configuration_label, group = configuration_label)
) +
  geom_hline(yintercept = 0, colour = "#6B747B", linewidth = 0.55) +
  geom_line(linewidth = 0.90) +
  geom_point(size = 3.0, stroke = 0.72, fill = paper) +
  scale_colour_manual(values = convergence_colours, name = NULL) +
  scale_shape_manual(values = convergence_shapes, name = NULL) +
  scale_linetype_manual(values = convergence_linetypes, name = NULL) +
  scale_y_continuous(
    breaks = seq(0, 1.0, by = 0.2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  coord_cartesian(ylim = c(-0.07, 1.06)) +
  labs(
    x = "Gaussian graph bandwidth (km)",
    y = "Estimate minus 465-site estimate (percentage points)"
  ) +
  jr_theme(10.0) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.25),
    legend.position = "none",
    legend.text = element_text(size = 8.0)
  )
save_figure("supp_spatial_convergence", convergence_figure, 8.3, 5.0)

###############################################################################
# Main-text Figure 6: held-out scale profile and annual recurrence.
###############################################################################

year_scale <- fread(file.path(
  confirmatory_dir, "confirmatory_year_scale_effects.csv"
))[day_definition == "utc" & grepl("^graph_h_", metric)]
metadata <- unique(fread(file.path(
  confirmatory_dir, "confirmatory_graph_metadata.csv"
))[definition_index == 1L, .(metric, bandwidth_km)])
year_scale <- merge(year_scale, metadata, by = "metric")
year_scale[, bandwidth_label := factor(
  sprintf("%s km", format(round(bandwidth_km), big.mark = ",")),
  levels = sprintf("%s km", format(round(sort(unique(bandwidth_km))),
                                   big.mark = ","))
)]
held_out_scale <- year_scale[analysis_role == "confirmatory"]
scale_summary <- fread(file.path(
  confirmatory_dir, "confirmatory_scale_results.csv"
))[day_definition == "utc" & grepl("^graph_h_", metric)]
scale_summary <- merge(scale_summary, metadata, by = "metric")
scale_summary[, bandwidth_label := factor(
  sprintf("%s km", format(round(bandwidth_km), big.mark = ",")),
  levels = levels(year_scale$bandwidth_label)
)]

figure3 <- ggplot() +
  annotate("rect", xmin = -40, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#F3F7FA", colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#6B747B") +
  geom_point(
    data = held_out_scale,
    aes(100 * yearly_relative_effect, bandwidth_label),
    position = position_jitter(width = 0, height = 0.13, seed = 20260809),
    shape = 21, size = 2.0, stroke = 0.25, colour = paper,
    fill = "#8FB6D1", alpha = 0.72
  ) +
  geom_errorbarh(
    data = scale_summary,
    aes(y = bandwidth_label, xmin = 100 * ci_lower, xmax = 100 * ci_upper),
    height = 0, linewidth = 1.05, colour = blue
  ) +
  geom_point(
    data = scale_summary,
    aes(100 * estimate, bandwidth_label),
    shape = 21, size = 4.0, stroke = 0.85, colour = blue, fill = paper
  ) +
  geom_text(
    data = scale_summary,
    aes(x = 11.5, y = bandwidth_label,
        label = sprintf("%d/33", negative_years)),
    hjust = 0.5, size = 3.2, colour = ink
  ) +
  scale_x_continuous(
    limits = c(-40, 16), breaks = seq(-40, 10, 10),
    labels = function(x) sprintf("%g%%", x)
  ) +
  labs(
    x = "Annual high-to-middle graph-dispersion contrast",
    y = "Gaussian graph bandwidth"
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.28),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )
save_figure("fig3_multiscale_evidence", figure3, 8.3, 4.8)

###############################################################################
# Main-text Figure 1: balanced power across unknown mechanisms.
###############################################################################

sim_summary <- fread(file.path(corrected_dir, "graph_simulation_summary.csv"))
method_columns <- c("rejection_graph", "rejection_variogram_profile",
                    "rejection_variance",
                    "rejection_nearest", "rejection_binned",
                    "rejection_moran", "rejection_geary")
method_labels <- c("Graph profile", "Five-bin variogram", "Spatial variance",
                   "Nearest-neighbour", "Single-bin semivariance",
                   "Moran's I", "Geary's C")
method_colours <- c(blue, "#B24783", orange, teal, gold, purple, "#65737E")

power <- melt(
  sim_summary[family == "alternative"],
  id.vars = c("scenario", "mechanism", "target_profile"),
  measure.vars = method_columns, variable.name = "method", value.name = "rate"
)
power[, mechanism := factor(mechanism,
                            levels = c("amplitude", "range", "gradient"),
                            labels = c("Amplitude", "Range", "Gradient"))]
power[, method := factor(method, levels = method_columns, labels = method_labels)]
power <- power[!grepl("_plus$", scenario)]
power[, strength := factor(
  fcase(
    grepl("_weak$", scenario), "Weak",
    grepl("_medium$", scenario), "Medium",
    grepl("_strong$", scenario), "Strong"
  ), levels = c("Weak", "Medium", "Strong")
)]
power_summary <- power[, .(
  minimum = min(rate), maximum = max(rate)
), by = method]
method_order <- power_summary[order(minimum)]$method
power[, method := factor(method, levels = method_order)]
power_summary[, method := factor(method, levels = method_order)]

figure4 <- ggplot() +
  geom_segment(
    data = power_summary,
    aes(x = 100 * minimum, xend = 100 * maximum, y = method, yend = method),
    colour = pale, linewidth = 3.6, lineend = "round"
  ) +
  geom_point(
    data = power,
    aes(100 * rate, method, fill = mechanism, shape = strength),
    size = 3.0, stroke = 0.45, colour = paper, alpha = 0.90
  ) +
  geom_point(
    data = power_summary,
    aes(100 * minimum, method),
    shape = 23, size = 3.5, stroke = 0.75, colour = ink, fill = paper
  ) +
  geom_text(
    data = power_summary,
    aes(100 * minimum, method, label = sprintf("%.1f", 100 * minimum)),
    nudge_x = -1.6, hjust = 1, size = 3.0, colour = ink
  ) +
  scale_fill_manual(values = c(Amplitude = blue, Range = teal,
                               Gradient = orange), name = "Mechanism") +
  scale_shape_manual(values = c(Weak = 21, Medium = 22, Strong = 24),
                     name = "Effect strength") +
  scale_x_continuous(
    limits = c(0, 102), breaks = seq(0, 100, 20),
    labels = function(x) sprintf("%g%%", x)
  ) +
  labs(
    x = "Rejection rate across nine alternatives",
    y = NULL
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.28),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )
save_figure("fig4_simulation_diagnostics", figure4, 8.3, 4.8)

###############################################################################
# Main-text Figure 2: coverage at the application record length.
###############################################################################

year_sim <- fread(file.path(corrected_dir,
                            "year_inference_simulation_summary.csv"))
coverage <- melt(
  year_sim[effect == "null" & innovation == "gaussian"],
  id.vars = c("sample_size", "rho"),
  measure.vars = c("coverage_student", "coverage_nw2", "coverage_hac"),
  variable.name = "method", value.name = "coverage"
)
coverage[, method := factor(
  method, levels = c("coverage_student", "coverage_nw2", "coverage_hac"),
  labels = c("Student", "Fixed lag 2", "Growing HAC")
)]
coverage[, row_label := sprintf("rho = %.1f    %s", rho, method)]
row_top_to_bottom <- unlist(lapply(
  c(0.6, 0.3, 0.0),
  function(r) sprintf("rho = %.1f    %s", r,
                      c("Student", "Fixed lag 2", "Growing HAC"))
))
coverage[, row_label := factor(row_label, levels = rev(row_top_to_bottom))]
coverage[, summer_label := factor(sample_size,
                                  levels = c(20, 33, 60, 120))]

figure4b <- ggplot(coverage,
                   aes(summer_label, row_label, fill = coverage)) +
  geom_tile(width = 0.92, height = 0.84, colour = paper, linewidth = 1.0) +
  geom_tile(
    data = coverage[sample_size == 33],
    width = 0.92, height = 0.84, fill = NA, colour = ink, linewidth = 0.9
  ) +
  geom_text(aes(label = sprintf("%.1f", 100 * coverage)),
            size = 3.35, colour = ink) +
  geom_hline(yintercept = c(3.5, 6.5), colour = grid, linewidth = 0.7) +
  scale_fill_gradientn(
    colours = c(orange, "#E9B58F", paper, light_blue, blue),
    values = scales::rescale(c(0.65, 0.78, 0.88, 0.94, 0.96)),
    limits = c(0.65, 0.96), oob = scales::squish,
    name = "Coverage"
  ) +
  labs(
    x = "Number of summers (outlined column is the 33-summer application)",
    y = NULL
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.key.width = grid::unit(28, "mm")
  )
save_figure("fig4b_asymptotic_diagnostics", figure4b, 8.3, 5.1)

###############################################################################
# Main-text Figure 7: robustness and estimand sensitivity.
###############################################################################
review_robustness <- fread(file.path(
  confirmatory_dir, "sensitivity_robustness_summary.csv"))
review_robustness[estimand == "log_effect", `:=`(
  estimate = transformed_effect,
  ci_lower = transformed_ci_lower,
  ci_upper = transformed_ci_upper
)]
review_robustness <- review_robustness[
  analysis %chin% c("equal_site_ratio", "area_weighted_fixed_labels",
                    "area_weighted_relabelled_primary_hours",
                    "exponential_matched_distance",
                    "compact_quadratic_matched_distance",
                    "site_month_climatology_anomaly",
                    "site_year_month_anomaly",
                    "site_year_month_linear_detrended",
                    "leave_one_year_daily_climatology_anomaly") |
    estimand == "log_effect",
  .(analysis = fcase(
      analysis == "equal_site_ratio" & estimand == "ratio_effect",
        "Primary ratio",
      analysis == "equal_site_ratio" & estimand == "log_effect",
        "Log ratio (back-transformed)",
      analysis == "area_weighted_fixed_labels", "Area-weighted graph",
      analysis == "area_weighted_relabelled_primary_hours",
        "Area-weighted graph; labels recomputed",
      analysis == "exponential_matched_distance", "Exponential kernel",
      analysis == "compact_quadratic_matched_distance", "Compact kernel",
      analysis == "site_month_climatology_anomaly", "Monthly-climatology anomaly",
      analysis == "site_year_month_anomaly", "Within-year monthly anomaly",
      analysis == "site_year_month_linear_detrended",
        "Site-record linear-detrended field",
      analysis == "leave_one_year_daily_climatology_anomaly",
        "Leave-one-year daily-climatology anomaly"),
    estimate, ci_lower, ci_upper,
    family = fcase(
      analysis == "equal_site_ratio", "Primary and transformation",
      grepl("weighted|matched", analysis), "Graph construction",
      default = "Field definition"))
]
dense <- fread(file.path(dense_dir, "dense_primary_results.csv"))[
  configuration %chin% c("dense_465_fixed_labels", "dense_465_recomputed"),
  .(analysis = fifelse(configuration == "dense_465_fixed_labels",
                       "465 sites; primary events",
                       "465 sites; events recomputed"),
    estimate, ci_lower, ci_upper)
]
robustness <- rbind(
  review_robustness,
  dense[, .(analysis, estimate, ci_lower, ci_upper,
            family = "Spatial resolution")], fill = TRUE
)
display_order <- c(
  "Primary ratio",
  "Log ratio (back-transformed)",
  "Area-weighted graph",
  "Area-weighted graph; labels recomputed",
  "Exponential kernel",
  "Compact kernel",
  "465 sites; primary events",
  "465 sites; events recomputed",
  "Site-record linear-detrended field",
  "Leave-one-year daily-climatology anomaly",
  "Monthly-climatology anomaly",
  "Within-year monthly anomaly"
)
robustness[, analysis := factor(analysis, levels = rev(display_order))]
robustness[, family := factor(
  family,
  levels = c("Primary and transformation", "Graph construction",
             "Spatial resolution", "Field definition")
)]

figure5 <- ggplot(robustness,
                  aes(100 * estimate, analysis, colour = family)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 4.5,
           fill = "#FBF5F0", colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#6B747B") +
  geom_errorbarh(aes(xmin = 100 * ci_lower, xmax = 100 * ci_upper),
                 height = 0, linewidth = 1.15) +
  geom_point(size = 3.7) +
  geom_point(
    data = robustness[analysis == "Primary ratio"],
    shape = 23, size = 5.0, stroke = 0.9, colour = ink, fill = paper
  ) +
  geom_text(
    aes(x = 4.2, label = sprintf("%+.1f", 100 * estimate)),
    hjust = 1, size = 3.1, colour = ink
  ) +
  geom_hline(yintercept = 4.5, colour = "#D6C5B8", linewidth = 0.7) +
  scale_colour_manual(values = c(
    "Primary and transformation" = ink, "Graph construction" = teal,
    "Field definition" = purple, "Spatial resolution" = orange),
    name = "Analysis family") +
  scale_x_continuous(
    limits = c(-16, 5), breaks = seq(-15, 5, 5),
    labels = function(x) sprintf("%g%%", x)
  ) +
  labs(
    x = "Mean profile effect (95% t-based interval)",
    y = NULL
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.28),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )
save_figure("fig5_application_robustness", figure5, 8.3, 5.8)

###############################################################################
# Main-text Figure 7: non-development measurement and effect agreement.
###############################################################################
measurement <- fread(file.path(
  noaa_extension_dir, "noaa_extension_era5_matched.csv.gz"
))
measurement_plot <- measurement[seq.int(1L, .N, by = 20L)]

measurement_panel <- ggplot(
  measurement_plot, aes(era_wbt_c, observed_wbt_c)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = "#6B747B", linewidth = 0.65) +
  geom_point(colour = blue, alpha = 0.075, size = 0.55) +
  annotate(
    "label", x = 1.0, y = 32.7, hjust = 0, vjust = 1,
    label = paste(
      "175,172 exact-hour matches",
      "Bias -0.26 degrees C; MAE 0.94; RMSE 1.27",
      sep = "\n"
    ),
    size = 3.0, colour = ink, fill = paper, label.size = 0.20
  ) +
  coord_equal(xlim = c(-1, 34), ylim = c(-1, 34), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = seq(0, 30, 10)) +
  labs(
    x = "ERA5-Land WBT (degrees C)",
    y = "NOAA station WBT (degrees C)",
    title = "A  Exact-hour measurement agreement"
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.28),
    panel.grid.major.y = element_line(colour = grid, linewidth = 0.28),
    legend.position = "none"
  )

station_graph <- fread(file.path(
  noaa_extension_dir, "noaa_extension_scale_summary.csv"
))
station_graph[, `:=`(
  observed_value = station_effect_estimate,
  era_value = era_effect_estimate,
  equality_value = 50 * (station_effect_estimate + era_effect_estimate),
  bandwidth_label = sprintf("%s km",
                            format(round(bandwidth_km), big.mark = ",")),
  label_x = fcase(
    bandwidth_km < 150, -10.4,
    bandwidth_km < 300, -15.2,
    bandwidth_km < 600, -17.0,
    bandwidth_km < 1500, -19.5,
    default = -25.2
  ),
  label_y = fcase(
    bandwidth_km < 150, -4.7,
    bandwidth_km < 300, -5.8,
    bandwidth_km < 600, -10.7,
    bandwidth_km < 1500, -15.7,
    default = -19.2
  )
)]

effect_panel <- ggplot(station_graph,
                       aes(100 * era_value, 100 * observed_value)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = "#6B747B", linewidth = 0.75) +
  geom_segment(
    aes(xend = equality_value, yend = equality_value),
    colour = pale, linewidth = 2.2, lineend = "round"
  ) +
  geom_point(
    aes(fill = bandwidth_km), shape = 21, size = 5.1,
    stroke = 0.85, colour = ink
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = bandwidth_label),
    size = 3.1, colour = ink, fontface = "bold"
  ) +
  annotate("text", x = -9, y = -6.5, label = "equality",
           angle = 45, hjust = 0, size = 3.0, colour = muted) +
  scale_fill_gradient(low = orange, high = blue, trans = "log10",
                      guide = "none") +
  scale_x_continuous(labels = function(x) sprintf("%g%%", x),
                     breaks = seq(-30, 10, 10)) +
  scale_y_continuous(labels = function(x) sprintf("%g%%", x),
                     breaks = seq(-30, 10, 10)) +
  coord_equal(xlim = c(-32, 8), ylim = c(-32, 8), expand = FALSE) +
  labs(
    x = "ERA5-Land high-to-middle graph contrast",
    y = "NOAA station high-to-middle graph contrast",
    title = "B  Effect agreement at common event times"
  ) +
  jr_theme(10.2) +
  theme(
    panel.grid.major.x = element_line(colour = grid, linewidth = 0.28),
    panel.grid.major.y = element_line(colour = grid, linewidth = 0.28),
    legend.position = "none"
  )

figure6 <- measurement_panel + effect_panel
save_figure("fig6_noaa_agreement", figure6, 8.3, 4.7)

# Copy each canonical vector file to the manuscript names used by LaTeX.  The
# smooth-map renderer below replaces Figures 3 and 4 after these copies.
portable_figures <- setNames(
  c(
    "figure01_simulation_diagnostics.pdf",
    "figure02_study_area.pdf",
    "figure05_energy_decomposition.pdf",
    "figure06_application_robustness.pdf",
    "figure07_noaa_agreement.pdf",
    "supp_dense_bandwidth_profile.pdf",
    "supp_spatial_convergence.pdf"
  ),
  c(
    file.path(output_dir, "fig4_simulation_diagnostics.pdf"),
    file.path(corrected_dir, "fig5_study_area.pdf"),
    file.path(output_dir, "fig_energy_decomposition.pdf"),
    file.path(output_dir, "fig5_application_robustness.pdf"),
    file.path(output_dir, "fig6_noaa_agreement.pdf"),
    file.path(output_dir, "supp_dense_bandwidth_profile.pdf"),
    file.path(output_dir, "supp_spatial_convergence.pdf")
  )
)
copy_status <- mapply(
  function(source, destination) {
    file.copy(
      source, file.path(manuscript_figure_dir, destination), overwrite = TRUE
    )
  },
  names(portable_figures), unname(portable_figures),
  USE.NAMES = FALSE
)
if (!all(copy_status)) stop("Failed to copy canonical publication figures")

# Finish with the primary-grid smooth map renderer.  It reuses `figure3` above
# as panel (D), adds the three exact 121-node scale-attribution maps as panels
# (A)-(C), and replaces the observed-node square map with a display-only
# thin-plate surface.  Keeping this call inside the main figure build prevents
# a later rerun of this script from restoring the superseded map versions.
primary_map_script <- file.path(code_dir, "53_primary_smooth_spatial_surfaces.R")
if (!file.exists(primary_map_script)) {
  stop("Missing primary smooth-map renderer: ", primary_map_script)
}
source(primary_map_script, local = TRUE)

cat(paste(
  "Wrote unified JRSS C main-text and supplementary figures as PDF files",
  "to output_jrssc.\n"
))
