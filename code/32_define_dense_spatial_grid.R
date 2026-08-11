###############################################################################
# Define a nested spatial-resolution sensitivity grid.
#
# The dense lattice halves both sampling intervals used by the frozen primary
# grid: 0.8 degrees longitude and 0.9 degrees latitude. Every original site is
# therefore retained exactly. This grid is a secondary resolution analysis and
# never replaces data/grid/eastern_china_121_sites.csv.
###############################################################################

library(data.table)
library(terra)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

archive <- file.path(project_dir, "data", "era5",
                     "era5_201506_eastern_china.zip")
primary_file <- file.path(project_dir, "data", "grid",
                          "eastern_china_121_sites.csv")
grid_dir <- file.path(project_dir, "data", "grid")
output_dir <- file.path(project_dir, "output_dense")
dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(archive)) stop("Missing cached reference archive: ", archive)
if (!file.exists(primary_file)) stop("Missing frozen primary manifest")

candidate <- CJ(
  lon = seq(105.0, 125.0, by = 0.8),
  lat = seq(42.0, 20.4, by = -0.9)
)
setorder(candidate, lon, lat)
candidate[, candidate_id := .I]
stopifnot(nrow(candidate) == 650L)

extract_dir <- tempfile("era5-dense-grid-audit-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
unzip(archive, files = "data.grib", exdir = extract_dir)
grib_path <- file.path(extract_dir, "data.grib")

x <- rast(grib_path)
temperature_index <- grep("2 metre temperature", names(x), fixed = TRUE)[1]
if (is.na(temperature_index)) stop("Could not find 2-m temperature")
points <- vect(candidate[, .(lon, lat)], geom = c("lon", "lat"),
               crs = "EPSG:4326")
reference_temperature <- extract(x[[temperature_index]], points,
                                 ID = FALSE)[[1]]
candidate[, era5_land_valid := is.finite(reference_temperature)]

coordinate_key <- function(lon, lat) sprintf("%.6f,%.6f", lon, lat)
primary <- fread(primary_file)
if (nrow(primary) != 121L) stop("Frozen primary manifest is not 121 sites")
primary[, coordinate_key := coordinate_key(lon, lat)]
candidate[, coordinate_key := coordinate_key(lon, lat)]
candidate[, original_site_id := primary$site_id[
  match(coordinate_key, primary$coordinate_key)
]]
candidate[, is_original_site := !is.na(original_site_id)]

if (candidate[, sum(is_original_site)] != 121L) {
  stop("Dense lattice does not contain all 121 primary coordinates")
}
if (candidate[is_original_site == TRUE, any(!era5_land_valid)]) {
  stop("A frozen primary location fails the dense-grid land mask")
}

manifest <- candidate[era5_land_valid == TRUE, .(
  lon, lat, original_site_id, is_original_site
)]
setorder(manifest, lon, lat)
manifest[, dense_site_id := .I]
setcolorder(manifest, c("dense_site_id", "lon", "lat",
                        "original_site_id", "is_original_site"))

if (manifest[, sum(is_original_site)] != 121L) {
  stop("Dense land manifest lost a primary location")
}
if (manifest[, anyDuplicated(coordinate_key(lon, lat))]) {
  stop("Duplicate coordinates in dense manifest")
}

fwrite(manifest, file.path(grid_dir, "eastern_china_dense_sites.csv"))
fwrite(candidate, file.path(output_dir, "dense_grid_audit.csv"))
summary <- data.table(
  longitude_spacing_degrees = 0.8,
  latitude_spacing_degrees = 0.9,
  candidate_sites = nrow(candidate),
  retained_land_sites = nrow(manifest),
  reused_primary_sites = manifest[, sum(is_original_site)],
  new_sites = manifest[, sum(!is_original_site)]
)
fwrite(summary, file.path(output_dir, "dense_grid_summary.csv"))

cat(sprintf("Dense candidate lattice: %d sites\n", nrow(candidate)))
cat(sprintf("ERA5-Land-valid dense sites: %d\n", nrow(manifest)))
cat(sprintf("Primary sites nested exactly: %d\n",
            manifest[, sum(is_original_site)]))
cat(sprintf("New CDS point requests required: %d\n",
            manifest[, sum(!is_original_site)]))
cat("Dense manifest: data/grid/eastern_china_dense_sites.csv\n")
