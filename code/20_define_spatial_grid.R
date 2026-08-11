###############################################################################
# Reconstruct and freeze the spatial sampling rule used by the corrected study.
#
# Rule: within 105--125 E and 20--42 N on the native 0.1-degree ERA5-Land
# grid, start at the north-west corner and retain every 16th longitude cell
# and every 18th latitude cell. This creates a 13 x 13 candidate lattice.
# A location is retained when the first cached 2-m-temperature field is not
# missing, which removes ERA5-Land ocean cells. The resulting 121-site
# manifest is compared with the historical 2015-06 coordinates and then frozen
# for every year and month.
###############################################################################

library(data.table)
library(terra)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

archive <- file.path(project_dir, "data", "era5",
                     "era5_201506_eastern_china.zip")
historical <- file.path(project_dir, "data", "era5",
                        "era5_201506_full.csv")
grid_dir <- file.path(project_dir, "data", "grid")
output_dir <- file.path(project_dir, "output_corrected")
dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(archive)) stop("Missing cached reference archive: ", archive)
if (!file.exists(historical)) stop("Missing historical coordinate file: ", historical)

candidate <- CJ(
  lon = seq(105.0, 124.2, by = 1.6),
  lat = seq(42.0, 20.4, by = -1.8)
)
setorder(candidate, lon, lat)
candidate[, candidate_id := .I]
stopifnot(nrow(candidate) == 169L)

extract_dir <- tempfile("era5-grid-audit-")
dir.create(extract_dir)
unzip(archive, files = "data.grib", exdir = extract_dir)
grib_path <- file.path(extract_dir, "data.grib")

x <- rast(grib_path)
temperature_index <- grep("2 metre temperature", names(x), fixed = TRUE)[1]
if (is.na(temperature_index)) stop("Could not find 2-m temperature in reference GRIB")

points <- vect(candidate[, .(lon, lat)], geom = c("lon", "lat"),
               crs = "EPSG:4326")
reference_temperature <- extract(x[[temperature_index]], points, ID = FALSE)[[1]]
candidate[, era5_land_valid := is.finite(reference_temperature)]

historical_sites <- unique(fread(historical, select = c("lon", "lat")))
coordinate_key <- function(lon, lat) sprintf("%.6f,%.6f", lon, lat)
historical_keys <- coordinate_key(historical_sites$lon, historical_sites$lat)
candidate[, historical_retained :=
            coordinate_key(lon, lat) %chin% historical_keys]

if (candidate[, any(era5_land_valid != historical_retained)]) {
  fwrite(candidate[era5_land_valid != historical_retained],
         file.path(output_dir, "spatial_grid_mismatches.csv"))
  stop("Reconstructed ERA5-Land mask does not reproduce historical sites")
}

manifest <- candidate[era5_land_valid == TRUE, .(lon, lat)]
setorder(manifest, lon, lat)
manifest[, site_id := .I]
setcolorder(manifest, c("site_id", "lon", "lat"))
stopifnot(nrow(manifest) == 121L)

fwrite(manifest, file.path(grid_dir, "eastern_china_121_sites.csv"))
fwrite(candidate, file.path(output_dir, "spatial_grid_audit.csv"))

cat(sprintf("Candidate lattice: %d sites\n", nrow(candidate)))
cat(sprintf("ERA5-Land-valid sites: %d\n", nrow(manifest)))
cat("Historical coordinate match: exact\n")
cat("Frozen manifest: data/grid/eastern_china_121_sites.csv\n")
