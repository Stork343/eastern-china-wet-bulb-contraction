###############################################################################
# Publication study-area and spatial-support figure.
###############################################################################

library(data.table)
library(ggplot2)
library(maps)
library(patchwork)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
output_dir <- file.path(project_dir, "output_corrected")

audit_file <- file.path(output_dir, "spatial_grid_audit.csv")
dense_file <- file.path(project_dir, "data", "grid",
                        "eastern_china_dense_sites.csv")
metadata_file <- file.path(project_dir, "output_confirmatory",
                           "confirmatory_graph_metadata.csv")
if (!all(file.exists(c(audit_file, dense_file, metadata_file)))) {
  stop("Spatial manifests or graph metadata are missing")
}
audit <- fread(audit_file)
dense <- fread(dense_file)
metadata <- unique(fread(metadata_file)[definition_index == 1L,
  .(h_factor, bandwidth_km)])
setorder(metadata, h_factor)
world <- map_data("world")
world <- world[world$long >= 102 & world$long <= 128 &
                 world$lat >= 17 & world$lat <= 45, ]

ink <- "#24303A"
blue <- "#1F5A85"
teal <- "#008C72"
gold <- "#E6A700"
paper <- "#FFFFFF"
grid <- "#D9DEE3"

map_theme <- theme_minimal(base_size = 9.2, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    plot.title = element_text(face = "bold", size = 9.8,
                              margin = margin(b = 4)),
    panel.grid.major = element_line(colour = grid, linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "#9AA3AA", fill = NA,
                                linewidth = 0.35),
    axis.title = element_text(size = 8.8),
    axis.text = element_text(size = 7.7),
    legend.position = "bottom",
    legend.title = element_text(size = 8.2),
    legend.text = element_text(size = 7.5),
    plot.margin = margin(4, 5, 4, 5)
  )

map_base <- function() {
  ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "#F1F3F4", colour = "#7C858C", linewidth = 0.25) +
    coord_quickmap(xlim = c(104.2, 125.4), ylim = c(19.3, 42.8),
                   expand = FALSE) +
    scale_x_continuous(breaks = c(105, 110, 115, 120, 125)) +
    scale_y_continuous(breaks = c(20, 25, 30, 35, 40)) +
    labs(x = "Longitude (degrees E)", y = "Latitude (degrees N)") +
    map_theme
}

p_primary <- map_base() +
  geom_point(data = audit[era5_land_valid == FALSE], aes(lon, lat),
             shape = 4, colour = "#AEB5BA", size = 1.2, stroke = 0.42) +
  geom_point(data = audit[era5_land_valid == TRUE], aes(lon, lat),
             shape = 21, fill = blue, colour = paper, size = 1.9,
             stroke = 0.25) +
  labs(title = "Primary 121-site lattice",
       subtitle = "169 candidates; ocean cells removed")

p_nested <- map_base() +
  geom_point(data = dense[is_original_site == FALSE], aes(lon, lat),
             shape = 21, fill = teal, colour = paper, size = 1.15,
             stroke = 0.12, alpha = 0.9) +
  geom_point(data = dense[is_original_site == TRUE], aes(lon, lat),
             shape = 21, fill = blue, colour = ink, size = 1.65,
             stroke = 0.28) +
  labs(title = "Nested 465-site lattice",
       subtitle = "Primary sites retained; 344 sites added")

primary <- audit[era5_land_valid == TRUE, .(lon, lat)]
reference <- primary[which.min((lon - 115)^2 + (lat - 31)^2)]
lat0 <- mean(primary$lat) * pi / 180
primary[, distance_km := sqrt(
  ((lon - reference$lon) * 111.32 * cos(lat0))^2 +
    ((lat - reference$lat) * 110.57)^2
)]

kernel_panel <- function(h, title) {
  d <- copy(primary)
  d[, weight := exp(-(distance_km^2) / (2 * h^2))]
  map_base() +
    geom_point(data = d, aes(lon, lat, fill = weight), shape = 21,
               colour = ink, size = 2.0, stroke = 0.18) +
    geom_point(data = reference, aes(lon, lat), shape = 23, fill = gold,
               colour = ink, size = 3.0, stroke = 0.45) +
    scale_fill_gradient(low = "#F3F5F6", high = blue, limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        name = "Gaussian edge weight") +
    labs(title = sprintf("%s (h = %.0f km)", title, h))
}

p_local <- kernel_panel(min(metadata$bandwidth_km), "Local graph support")
p_broad <- kernel_panel(max(metadata$bandwidth_km), "Domain-scale graph support")

figure <- ((p_primary | p_nested) / (p_local | p_broad)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.tag = element_text(face = "bold",
                                                        colour = ink))) &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig5_study_area.pdf"), figure,
       width = 8.3, height = 7.1, device = cairo_pdf, bg = paper)
cat("Study-area and spatial-support figure written to output_corrected.\n")
