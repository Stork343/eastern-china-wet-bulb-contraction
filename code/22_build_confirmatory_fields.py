#!/usr/bin/env python3
"""Quality-control ERA5-Land point panels and build synchronous daily fields."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "data" / "era5_confirmatory" / "hourly_points"
OUTPUT_DIR = PROJECT_DIR / "data" / "era5_confirmatory" / "daily_fields"
AUDIT_FILE = PROJECT_DIR / "output_corrected" / "confirmatory_data_audit.csv"
DISCOVERY_YEARS = {2015, 2022}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--years", type=int, nargs="+")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def saturation_vapor_pressure_pa(temperature_k: np.ndarray) -> np.ndarray:
    temperature_c = temperature_k - 273.15
    return 611.2 * np.exp(17.67 * temperature_c / (temperature_c + 243.5))


def bolton_theta_e(temperature_k: np.ndarray, dewpoint_k: np.ndarray,
                   pressure_pa: np.ndarray) -> np.ndarray:
    vapor_pressure = saturation_vapor_pressure_pa(dewpoint_k)
    mixing_ratio = 0.622 * vapor_pressure / (pressure_pa - vapor_pressure)
    t_lcl = 1 / (1 / (dewpoint_k - 56) +
                 np.log(temperature_k / dewpoint_k) / 800) + 56
    theta_l = (temperature_k *
               (100000 / (pressure_pa - vapor_pressure)) ** 0.2854 *
               (temperature_k / t_lcl) ** (0.28 * mixing_ratio))
    return theta_l * np.exp((3036 / t_lcl - 1.78) * mixing_ratio *
                            (1 + 0.448 * mixing_ratio))


def bolton_wbt_c(temperature_k: np.ndarray, dewpoint_k: np.ndarray,
                 pressure_pa: np.ndarray, iterations: int = 40) -> np.ndarray:
    dewpoint_k = np.minimum(dewpoint_k, temperature_k)
    target = bolton_theta_e(temperature_k, dewpoint_k, pressure_pa)
    lower = dewpoint_k.copy()
    upper = temperature_k.copy()
    for _ in range(iterations):
        midpoint = (lower + upper) / 2
        saturated = bolton_theta_e(midpoint, midpoint, pressure_pa)
        move_lower = saturated < target
        lower = np.where(move_lower, midpoint, lower)
        upper = np.where(move_lower, upper, midpoint)
    return (lower + upper) / 2 - 273.15


def stull_wbt_c(temperature_k: np.ndarray,
                dewpoint_k: np.ndarray) -> np.ndarray:
    """Stull (2011) wet-bulb approximation used only as a sensitivity."""
    temperature_c = temperature_k - 273.15
    relative_humidity = 100 * (
        saturation_vapor_pressure_pa(dewpoint_k) /
        saturation_vapor_pressure_pa(temperature_k)
    )
    relative_humidity = np.clip(relative_humidity, 0, 100)
    return (
        temperature_c * np.arctan(
            0.151977 * np.sqrt(relative_humidity + 8.313659)
        ) +
        np.arctan(temperature_c + relative_humidity) -
        np.arctan(relative_humidity - 1.676331) +
        0.00391838 * relative_humidity ** 1.5 *
        np.arctan(0.023101 * relative_humidity) - 4.686035
    )


def validate_panel(frame: pd.DataFrame, year: int) -> dict[str, float | int | str]:
    required = {"time", "site_id", "requested_lon", "requested_lat",
                "t2m", "d2m", "sp"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{year}: missing columns {sorted(missing)}")
    if frame["time"].nunique() != 2232 or frame["site_id"].nunique() != 121:
        raise ValueError(f"{year}: incomplete time-site panel")
    counts = frame.groupby("time", observed=True)["site_id"].nunique()
    if not (counts == 121).all():
        raise ValueError(f"{year}: incomplete hourly spatial field")
    for variable in ("t2m", "d2m", "sp"):
        if not np.isfinite(frame[variable]).all():
            raise ValueError(f"{year}: nonfinite values in {variable}")
    if not frame["t2m"].between(220, 340).all():
        raise ValueError(f"{year}: implausible 2-m temperature")
    if not frame["sp"].between(45000, 110000).all():
        raise ValueError(f"{year}: implausible surface pressure")
    return {
        "year": year,
        "analysis_role": "discovery" if year in DISCOVERY_YEARS else "confirmatory",
        "hours_with_buffer": frame["time"].nunique(),
        "sites": frame["site_id"].nunique(),
        "dewpoint_exceedances_clipped": int((frame["d2m"] > frame["t2m"]).sum()),
        "t2m_min_c": float(frame["t2m"].min() - 273.15),
        "t2m_max_c": float(frame["t2m"].max() - 273.15),
        "sp_min_pa": float(frame["sp"].min()),
        "sp_max_pa": float(frame["sp"].max()),
    }


def select_peak_fields(frame: pd.DataFrame, year: int,
                       definition: str, wbt_column: str = "wbt") -> pd.DataFrame:
    working = frame.copy()
    if definition in {"utc", "utc_stull", "sitewise_max"}:
        working["analysis_date"] = working["time"].dt.floor("D")
        in_jja = working["analysis_date"].dt.month.isin((6, 7, 8))
    elif definition == "utc_plus_8":
        local_time = working["time"] + pd.Timedelta(hours=8)
        working["analysis_date"] = local_time.dt.floor("D")
        in_jja = working["analysis_date"].dt.month.isin((6, 7, 8))
    else:
        raise ValueError(definition)
    working = working.loc[in_jja].copy()
    working["analysis_wbt"] = working[wbt_column]

    if definition == "sitewise_max":
        selected = (working.sort_values(
            ["analysis_date", "site_id", "analysis_wbt", "time"],
            ascending=[True, True, False, True]
        ).drop_duplicates(["analysis_date", "site_id"]))
        regional = (selected.groupby("analysis_date", observed=True)
                    ["analysis_wbt"].mean().rename("regional_mean_wbt"))
        selected = selected.join(regional, on="analysis_date")
    else:
        regional = (working.groupby("time", observed=True)["analysis_wbt"]
                    .mean().rename("regional_mean_wbt").reset_index())
        regional["analysis_date"] = (
            regional["time"].dt.floor("D") if definition != "utc_plus_8"
            else (regional["time"] + pd.Timedelta(hours=8)).dt.floor("D")
        )
        peak_rows = (regional.sort_values(
            ["analysis_date", "regional_mean_wbt", "time"],
            ascending=[True, False, True]
        ).drop_duplicates("analysis_date"))
        selected = working.merge(
            peak_rows[["analysis_date", "time", "regional_mean_wbt"]],
            on=["analysis_date", "time"], how="inner", validate="many_to_one"
        )

    if selected["analysis_date"].nunique() != 92 or len(selected) != 92 * 121:
        raise ValueError(f"{year}: incomplete {definition} daily fields")
    selected["wbt"] = selected["analysis_wbt"]
    selected["wbt_method"] = "Stull_2011" if wbt_column == "wbt_stull" else \
        "pressure_aware_Bolton"
    selected["day_definition"] = definition
    selected["year"] = year
    selected["month"] = selected["analysis_date"].dt.month
    selected["record_id"] = selected["analysis_date"].dt.strftime("%Y%m")
    selected["analysis_role"] = (
        "discovery" if year in DISCOVERY_YEARS else "confirmatory"
    )
    return selected


def process_year(path: Path, overwrite: bool) -> tuple[Path, dict]:
    import xarray as xr

    year = int(path.name.split("_")[2])
    output = OUTPUT_DIR / f"era5_land_{year}_jja_daily_fields.csv.gz"
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}")
        frame = xr.open_dataset(path, engine="h5netcdf").to_dataframe().reset_index()
        frame["time"] = pd.to_datetime(frame["time"], utc=True).dt.tz_localize(None)
        return output, validate_panel(frame, year)

    dataset = xr.open_dataset(path, engine="h5netcdf")
    frame = dataset.to_dataframe().reset_index()
    frame["time"] = pd.to_datetime(frame["time"], utc=True).dt.tz_localize(None)
    audit = validate_panel(frame, year)
    frame["d2m"] = np.minimum(frame["d2m"], frame["t2m"])
    frame["wbt"] = bolton_wbt_c(
        frame["t2m"].to_numpy(), frame["d2m"].to_numpy(),
        frame["sp"].to_numpy()
    )
    frame["wbt_stull"] = stull_wbt_c(
        frame["t2m"].to_numpy(), frame["d2m"].to_numpy()
    )
    if not np.isfinite(frame["wbt"]).all() or not frame["wbt"].between(-60, 50).all():
        raise ValueError(f"{year}: invalid wet-bulb temperatures")
    if not np.isfinite(frame["wbt_stull"]).all():
        raise ValueError(f"{year}: invalid Stull wet-bulb temperatures")

    fields = pd.concat([
        select_peak_fields(frame, year, "utc"),
        select_peak_fields(frame, year, "utc_plus_8"),
        select_peak_fields(frame, year, "utc_stull", "wbt_stull"),
        select_peak_fields(frame, year, "sitewise_max"),
    ], ignore_index=True)
    keep = [
        "analysis_role", "year", "month", "record_id", "analysis_date",
        "day_definition", "wbt_method", "time", "site_id", "requested_lon", "requested_lat",
        "regional_mean_wbt", "wbt", "t2m", "d2m", "sp",
    ]
    keep.extend(variable for variable in ("u10", "v10", "swvl1", "ssrd")
                if variable in fields.columns)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fields[keep].to_csv(output, index=False, compression="gzip")
    audit["wbt_min_c"] = float(fields["wbt"].min())
    audit["wbt_max_c"] = float(fields["wbt"].max())
    print(f"[DONE] {output.name}: 4 definitions x 92 days x 121 sites")
    return output, audit


def main() -> None:
    args = parse_args()
    paths = sorted(INPUT_DIR.glob("era5_land_*_jja_121sites.nc"))
    if args.years:
        requested = set(args.years)
        paths = [path for path in paths if int(path.name.split("_")[2]) in requested]
    if not paths:
        raise SystemExit(f"No yearly point panels found in {INPUT_DIR}")
    audits = [process_year(path, args.overwrite)[1] for path in paths]
    AUDIT_FILE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audits).sort_values("year").to_csv(AUDIT_FILE, index=False)
    print(f"Audit: {AUDIT_FILE}")


if __name__ == "__main__":
    main()
