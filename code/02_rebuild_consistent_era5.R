###############################################################################
# 02_rebuild_consistent_era5.R
#
# Rebuild the six ERA5-Land monthly CSV files on one fixed spatial grid.
# The historical CSVs were produced with two different sampling schemes, so
# site_id was not a stable spatial key. This script preserves the original
# files and writes corrected data to data/era5_consistent/.
###############################################################################

library(data.table)
library(terra)

set.seed(20260802)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
source_dir <- file.path(project_dir, "data", "era5")
output_dir <- file.path(project_dir, "data", "era5_consistent")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

months <- c("201506", "201507", "201508", "202206", "202207", "202208")

# Script 20 reconstructs the 13 x 13 candidate lattice and ERA5-Land ocean
# omission rule from the cached GRIB, then freezes the 121 retained sites.
reference_file <- file.path(project_dir, "data", "grid",
                            "eastern_china_121_sites.csv")
if (!file.exists(reference_file)) {
  stop("Run 20_define_spatial_grid.R before rebuilding the monthly files")
}
reference_sites <- fread(reference_file)
setorder(reference_sites, site_id)
stopifnot(nrow(reference_sites) == 121L)

stull_wbt_c <- function(t_k, td_k) {
  t_c <- t_k - 273.15
  td_c <- td_k - 273.15
  e <- 6.112 * exp(17.67 * td_c / (td_c + 243.5))
  es <- 6.112 * exp(17.67 * t_c / (t_c + 243.5))
  rh <- pmin(pmax(100 * e / es, 1), 100)

  t_c * atan(0.151977 * sqrt(rh + 8.313659)) +
    atan(t_c + rh) - atan(rh - 1.676331) +
    0.00391838 * rh^(3 / 2) * atan(0.023101 * rh) - 4.686035
}

saturation_vapor_pressure_pa <- function(t_k) {
  t_c <- t_k - 273.15
  611.2 * exp(17.67 * t_c / (t_c + 243.5))
}

bolton_theta_e <- function(t_k, td_k, pressure_pa) {
  vapor_pressure <- saturation_vapor_pressure_pa(td_k)
  mixing_ratio <- 0.622 * vapor_pressure / (pressure_pa - vapor_pressure)
  t_lcl <- 1 / (1 / (td_k - 56) + log(t_k / td_k) / 800) + 56
  theta_l <- t_k * (100000 / (pressure_pa - vapor_pressure))^0.2854 *
    (t_k / t_lcl)^(0.28 * mixing_ratio)
  theta_l * exp((3036 / t_lcl - 1.78) * mixing_ratio *
                  (1 + 0.448 * mixing_ratio))
}

# Pressure-aware thermodynamic wet-bulb approximation. We calculate the
# parcel's Bolton equivalent potential temperature, then solve for the
# saturated temperature at the original surface pressure (Normand's rule).
bolton_wbt_c <- function(t_k, td_k, pressure_pa, iterations = 40L) {
  td_k <- pmin(td_k, t_k)
  target <- bolton_theta_e(t_k, td_k, pressure_pa)
  lower <- td_k
  upper <- t_k
  for (i in seq_len(iterations)) {
    midpoint <- (lower + upper) / 2
    saturated_theta_e <- bolton_theta_e(midpoint, midpoint, pressure_pa)
    move_lower <- saturated_theta_e < target
    lower[move_lower] <- midpoint[move_lower]
    upper[!move_lower] <- midpoint[!move_lower]
  }
  (lower + upper) / 2 - 273.15
}

select_layers <- function(layer_names, pattern, expected) {
  idx <- grep(pattern, layer_names, fixed = TRUE)
  if (length(idx) != expected) {
    stop(sprintf("Expected %d layers matching '%s'; found %d",
                 expected, pattern, length(idx)))
  }
  idx
}

rebuild_month <- function(ym) {
  zip_path <- file.path(source_dir, sprintf("era5_%s_eastern_china.zip", ym))
  if (!file.exists(zip_path)) stop("Missing cached archive: ", zip_path)

  month_tmp <- tempfile(sprintf("era5-%s-", ym))
  dir.create(month_tmp)
  on.exit(unlink(month_tmp, recursive = TRUE, force = TRUE), add = TRUE)
  unzip(zip_path, files = "data.grib", exdir = month_tmp)

  grib_path <- file.path(month_tmp, "data.grib")
  x <- rast(grib_path)
  layer_names <- names(x)
  n_time <- nlyr(x) %/% 6L
  if (n_time * 6L != nlyr(x)) stop("Unexpected number of GRIB layers: ", nlyr(x))

  idx_t2m <- select_layers(layer_names, "2 metre temperature", n_time)
  idx_d2m <- select_layers(layer_names, "2 metre dewpoint temperature", n_time)
  idx_sp <- select_layers(layer_names, "Surface pressure", n_time)
  idx_u10 <- select_layers(layer_names, "10 metre u wind component", n_time)
  idx_v10 <- select_layers(layer_names, "10 metre v wind component", n_time)
  idx_ssrd <- select_layers(layer_names, "Surface solar radiation downwards", n_time)

  points <- vect(reference_sites[, .(lon, lat)],
                 geom = c("lon", "lat"), crs = "EPSG:4326")
  vals <- extract(x, points, ID = FALSE)

  flatten <- function(idx) as.vector(as.matrix(vals[, idx, drop = FALSE]))
  t2m_k <- flatten(idx_t2m)
  # Interpolation occasionally places dew point up to 0.002 C above air
  # temperature. Enforce the physical saturation bound before all humidity
  # calculations and in the exported column.
  d2m_k <- pmin(flatten(idx_d2m), t2m_k)
  sp <- flatten(idx_sp)
  u10 <- flatten(idx_u10)
  v10 <- flatten(idx_v10)

  layer_time <- time(x)[idx_t2m]
  if (anyNA(layer_time)) stop("Missing GRIB timestamps for ", ym)

  out <- data.table(
    site_id = rep(reference_sites$site_id, times = n_time),
    lon = rep(reference_sites$lon, times = n_time),
    lat = rep(reference_sites$lat, times = n_time),
    time = rep(format(layer_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
               each = nrow(reference_sites)),
    wbt = bolton_wbt_c(t2m_k, d2m_k, sp),
    wbt_stull = stull_wbt_c(t2m_k, d2m_k),
    t2m = t2m_k - 273.15,
    d2m = d2m_k - 273.15,
    sp = sp,
    u10 = u10,
    v10 = v10,
    wspd = sqrt(u10^2 + v10^2),
    ssrd = flatten(idx_ssrd)
  )

  if (out[, uniqueN(paste(lon, lat))] != nrow(reference_sites)) {
    stop("Coordinate loss while rebuilding ", ym)
  }
  if (out[, uniqueN(time)] != n_time) stop("Time loss while rebuilding ", ym)
  if (out[, anyNA(wbt)]) stop("Missing WBT values while rebuilding ", ym)

  output_file <- file.path(output_dir, sprintf("era5_%s_consistent.csv", ym))
  fwrite(out, output_file)
  cat(sprintf("%s: %d sites x %d hours -> %s\n",
              ym, nrow(reference_sites), n_time, output_file))
  invisible(output_file)
}

for (ym in months) rebuild_month(ym)

# Cross-file invariants required by the spatial analysis.
rebuilt_files <- file.path(output_dir, sprintf("era5_%s_consistent.csv", months))
site_sets <- lapply(rebuilt_files, function(f) {
  unique(fread(f, select = c("site_id", "lon", "lat")))
})
reference_key <- site_sets[[1]][order(site_id)]
for (i in seq_along(site_sets)) {
  stopifnot(identical(reference_key, site_sets[[i]][order(site_id)]))
}

cat("All six files share an identical 121-site coordinate grid.\n")
