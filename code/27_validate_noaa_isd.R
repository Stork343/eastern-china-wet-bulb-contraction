###############################################################################
# Independent validation against NOAA Integrated Surface Database observations.
###############################################################################

library(data.table)
library(terra)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
era5_dir <- file.path(project_dir, "data", "era5")
isd_dir <- file.path(project_dir, "data", "noaa_isd_validation")
output_dir <- file.path(project_dir, "output_corrected")

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

bolton_wbt_c <- function(t_k, td_k, pressure_pa, iterations = 40L) {
  result <- rep(NA_real_, length(t_k))
  valid <- is.finite(t_k) & is.finite(td_k) & is.finite(pressure_pa) &
    pressure_pa > 45000 & pressure_pa < 110000
  if (!any(valid)) return(result)
  t_valid <- t_k[valid]
  td_valid <- pmin(td_k[valid], t_valid)
  pressure_valid <- pressure_pa[valid]
  target <- bolton_theta_e(t_valid, td_valid, pressure_valid)
  lower <- td_valid
  upper <- t_valid
  for (i in seq_len(iterations)) {
    midpoint <- (lower + upper) / 2
    saturated <- bolton_theta_e(midpoint, midpoint, pressure_valid)
    move_lower <- saturated < target
    lower[move_lower] <- midpoint[move_lower]
    upper[!move_lower] <- midpoint[!move_lower]
  }
  result[valid] <- (lower + upper) / 2 - 273.15
  result
}

select_layers <- function(layer_names, pattern, expected) {
  index <- grep(pattern, layer_names, fixed = TRUE)
  if (length(index) != expected) {
    stop("Expected ", expected, " layers matching ", pattern,
         "; found ", length(index))
  }
  index
}

observation_files <- file.path(isd_dir,
                               sprintf("noaa_isd_%d_jja.csv.gz", c(2015, 2022)))
if (any(!file.exists(observation_files))) {
  stop("Run 26_download_noaa_isd.py before validation")
}
observations <- rbindlist(lapply(observation_files, fread), fill = TRUE)
observations[, `:=`(
  STATION = as.character(STATION),
  time_utc = as.POSIXct(DATE, tz = "UTC"),
  year = as.integer(substr(DATE, 1, 4))
)]
accepted_qc <- c("0", "1", "4", "5", "9", "A", "C", "I", "M", "P", "R", "U")
observations <- observations[
  is.finite(temperature_c) & is.finite(dewpoint_c) &
    is.finite(station_pressure_pa) &
    as.character(temperature_qc) %chin% accepted_qc &
    as.character(dewpoint_qc) %chin% accepted_qc &
    as.character(slp_qc) %chin% accepted_qc
]
observations[, dewpoint_c := pmin(dewpoint_c, temperature_c)]
# Multiple source reports can occur at the same nominal station hour. Retain
# one deterministic record so high-frequency airports cannot duplicate hours.
setorder(observations, STATION, time_utc, REPORT_TYPE)
observations <- unique(observations, by = c("STATION", "time_utc"))
stations <- unique(observations[, .(
  STATION, LATITUDE, LONGITUDE, ELEVATION, NAME
)])

extract_month <- function(year, month) {
  ym <- sprintf("%d%02d", year, month)
  cache <- file.path(isd_dir, sprintf("era5_at_isd_stations_%s.csv.gz", ym))
  if (file.exists(cache)) return(fread(cache))
  archive <- file.path(era5_dir, sprintf("era5_%s_eastern_china.zip", ym))
  if (!file.exists(archive)) stop("Missing ERA5 archive: ", archive)
  extraction_dir <- tempfile(paste0("era5-isd-", ym, "-"))
  dir.create(extraction_dir)
  unzip(archive, files = "data.grib", exdir = extraction_dir)
  x <- rast(file.path(extraction_dir, "data.grib"))
  n_time <- nlyr(x) %/% 6L
  t_index <- select_layers(names(x), "2 metre temperature", n_time)
  d_index <- select_layers(names(x), "2 metre dewpoint temperature", n_time)
  p_index <- select_layers(names(x), "Surface pressure", n_time)

  points <- vect(stations[, .(lon = LONGITUDE, lat = LATITUDE)],
                 geom = c("lon", "lat"), crs = "EPSG:4326")
  values <- extract(x, points, ID = FALSE)
  flatten <- function(index) as.vector(as.matrix(values[, index, drop = FALSE]))
  t_k <- flatten(t_index)
  d_k <- pmin(flatten(d_index), t_k)
  pressure <- flatten(p_index)
  timestamps <- time(x)[t_index]
  if (anyNA(timestamps)) stop("Missing GRIB timestamps in ", ym)

  result <- data.table(
    STATION = rep(stations$STATION, times = n_time),
    time_utc = rep(as.POSIXct(timestamps, tz = "UTC"), each = nrow(stations)),
    era5_temperature_c = t_k - 273.15,
    era5_dewpoint_c = d_k - 273.15,
    era5_surface_pressure_pa = pressure,
    era5_wbt_c = bolton_wbt_c(t_k, d_k, pressure)
  )
  fwrite(result, cache)
  rm(x, values)
  unlink(extraction_dir, recursive = TRUE)
  result
}

era5 <- rbindlist(lapply(c(2015L, 2022L), function(year) {
  rbindlist(lapply(6:8, function(month) extract_month(year, month)))
}))
era5[, STATION := as.character(STATION)]

matched <- merge(observations, era5, by = c("STATION", "time_utc"))
matched <- matched[
  is.finite(era5_temperature_c) & is.finite(era5_dewpoint_c) &
    is.finite(era5_surface_pressure_pa) & is.finite(era5_wbt_c)
]
matched[, `:=`(
  observed_wbt_c = bolton_wbt_c(temperature_c + 273.15,
                                 dewpoint_c + 273.15,
                                 station_pressure_pa),
  observed_wbt_era5_pressure_c = bolton_wbt_c(
    temperature_c + 273.15, dewpoint_c + 273.15,
    era5_surface_pressure_pa
  )
)]

metric_summary <- function(observed, modeled) {
  complete <- is.finite(observed) & is.finite(modeled)
  if (sum(complete) < 2L) {
    return(list(n = sum(complete), bias = NA_real_, mae = NA_real_,
                rmse = NA_real_, correlation = NA_real_))
  }
  observed <- observed[complete]
  modeled <- modeled[complete]
  error <- modeled - observed
  list(
    n = length(observed),
    bias = mean(error, na.rm = TRUE),
    mae = mean(abs(error), na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    correlation = cor(observed, modeled, use = "complete.obs")
  )
}

station_results <- matched[, c(
  list(station_name = first(NAME), latitude = first(LATITUDE),
       longitude = first(LONGITUDE), elevation_m = first(ELEVATION)),
  setNames(metric_summary(observed_wbt_c, era5_wbt_c),
           paste0("wbt_", names(metric_summary(observed_wbt_c, era5_wbt_c))))
), by = .(year, STATION)]

overall <- rbindlist(list(
  matched[, c(list(variable = "temperature"),
              metric_summary(temperature_c, era5_temperature_c))],
  matched[, c(list(variable = "dewpoint"),
              metric_summary(dewpoint_c, era5_dewpoint_c))],
  matched[, c(list(variable = "wet_bulb_station_pressure"),
              metric_summary(observed_wbt_c, era5_wbt_c))],
  matched[, c(list(variable = "wet_bulb_common_era5_pressure"),
              metric_summary(observed_wbt_era5_pressure_c, era5_wbt_c))]
))

overall_by_year <- rbindlist(lapply(sort(unique(matched$year)), function(y) {
  z <- matched[year == y]
  rbindlist(list(
    z[, c(list(year = y, variable = "temperature"),
          metric_summary(temperature_c, era5_temperature_c))],
    z[, c(list(year = y, variable = "dewpoint"),
          metric_summary(dewpoint_c, era5_dewpoint_c))],
    z[, c(list(year = y, variable = "wet_bulb_station_pressure"),
          metric_summary(observed_wbt_c, era5_wbt_c))],
    z[, c(list(year = y, variable = "wet_bulb_common_era5_pressure"),
          metric_summary(observed_wbt_era5_pressure_c, era5_wbt_c))]
  ))
}))

fwrite(matched, file.path(output_dir, "noaa_isd_era5_matched.csv.gz"))
fwrite(station_results, file.path(output_dir, "noaa_isd_station_validation.csv"))
fwrite(overall, file.path(output_dir, "noaa_isd_overall_validation.csv"))
fwrite(overall_by_year,
       file.path(output_dir, "noaa_isd_year_validation.csv"))
cat(sprintf("Matched %d station-hour pairs across %d stations.\n",
            nrow(matched), uniqueN(matched$STATION)))
print(overall, digits = 4)
