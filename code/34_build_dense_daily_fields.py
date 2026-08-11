#!/usr/bin/env python3
"""Build daily fields for the 465-site nested-grid sensitivity analysis.

Two definitions are retained. ``primary_grid_peak`` uses the UTC peak times
and label-generating regional means from the frozen 121-site fields.
``dense_grid_peak`` selects the peak time and labels from the 465-site grid.
This separates spatial-support sensitivity from a full grid reanalysis.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
CODE_DIR = Path(__file__).resolve().parent
INPUT_DIR = PROJECT_DIR / "data" / "era5_dense" / "hourly_points"
OUTPUT_DIR = PROJECT_DIR / "data" / "era5_dense" / "daily_fields"
PRIMARY_DAILY_DIR = PROJECT_DIR / "data" / "era5_confirmatory" / "daily_fields"
AUDIT_FILE = PROJECT_DIR / "output_dense" / "dense_data_audit.csv"
DISCOVERY_YEARS = {2015, 2022}
EXPECTED_SITES = 465


def load_primary_builder():
    path = CODE_DIR / "22_build_confirmatory_fields.py"
    spec = importlib.util.spec_from_file_location("primary_field_builder", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRIMARY = load_primary_builder()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--years", type=int, nargs="+")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def validate_panel(frame: pd.DataFrame, year: int) -> dict[str, float | int | str]:
    required = {
        "time", "site_id", "requested_lon", "requested_lat", "t2m", "d2m", "sp"
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{year}: missing columns {sorted(missing)}")
    if frame.time.nunique() != 2232 or frame.site_id.nunique() != EXPECTED_SITES:
        raise ValueError(f"{year}: incomplete dense time-site panel")
    counts = frame.groupby("time", observed=True).site_id.nunique()
    if not (counts == EXPECTED_SITES).all():
        raise ValueError(f"{year}: incomplete dense hourly field")
    for variable in ("t2m", "d2m", "sp"):
        if not np.isfinite(frame[variable]).all():
            raise ValueError(f"{year}: nonfinite values in {variable}")
    if not frame.t2m.between(220, 340).all():
        raise ValueError(f"{year}: implausible 2-m temperature")
    if not frame.sp.between(45000, 110000).all():
        raise ValueError(f"{year}: implausible surface pressure")
    return {
        "year": year,
        "analysis_role": "discovery" if year in DISCOVERY_YEARS else "confirmatory",
        "hours_with_buffer": int(frame.time.nunique()),
        "sites": int(frame.site_id.nunique()),
        "dewpoint_exceedances_clipped": int((frame.d2m > frame.t2m).sum()),
        "t2m_min_c": float(frame.t2m.min() - 273.15),
        "t2m_max_c": float(frame.t2m.max() - 273.15),
        "sp_min_pa": float(frame.sp.min()),
        "sp_max_pa": float(frame.sp.max()),
    }


def utc_jja(frame: pd.DataFrame) -> pd.DataFrame:
    working = frame.copy()
    working["analysis_date"] = working.time.dt.floor("D")
    return working.loc[working.analysis_date.dt.month.isin((6, 7, 8))].copy()


def finalize(selected: pd.DataFrame, year: int, definition: str,
             label_mean: pd.Series | np.ndarray) -> pd.DataFrame:
    selected = selected.copy()
    selected["label_mean_wbt"] = np.asarray(label_mean)
    selected["wbt_method"] = "pressure_aware_Bolton"
    selected["analysis_definition"] = definition
    selected["year"] = year
    selected["month"] = selected.analysis_date.dt.month
    selected["record_id"] = selected.analysis_date.dt.strftime("%Y%m")
    selected["analysis_role"] = (
        "discovery" if year in DISCOVERY_YEARS else "confirmatory"
    )
    return selected


def select_dense_peak(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    working = utc_jja(frame)
    regional = (
        working.groupby("time", observed=True).wbt.mean()
        .rename("regional_mean_wbt").reset_index()
    )
    regional["analysis_date"] = regional.time.dt.floor("D")
    peaks = (
        regional.sort_values(
            ["analysis_date", "regional_mean_wbt", "time"],
            ascending=[True, False, True],
        ).drop_duplicates("analysis_date")
    )
    selected = working.merge(
        peaks[["analysis_date", "time", "regional_mean_wbt"]],
        on=["analysis_date", "time"], how="inner", validate="many_to_one",
    )
    if selected.analysis_date.nunique() != 92 or \
            len(selected) != 92 * EXPECTED_SITES:
        raise ValueError(f"{year}: incomplete dense-grid peak fields")
    return finalize(
        selected, year, "dense_grid_peak", selected.regional_mean_wbt
    )


def primary_peak_schedule(year: int) -> pd.DataFrame:
    path = PRIMARY_DAILY_DIR / f"era5_land_{year}_jja_daily_fields.csv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing primary daily field: {path}")
    schedule = pd.read_csv(
        path,
        usecols=["analysis_date", "day_definition", "time", "regional_mean_wbt"],
        parse_dates=["analysis_date", "time"],
    )
    schedule = schedule.loc[schedule.day_definition == "utc"].drop_duplicates(
        ["analysis_date", "time"]
    )
    if len(schedule) != 92 or schedule.analysis_date.nunique() != 92:
        raise ValueError(f"{year}: invalid primary peak schedule")
    return schedule.rename(columns={
        "regional_mean_wbt": "primary_regional_mean_wbt"
    })[["analysis_date", "time", "primary_regional_mean_wbt"]]


def select_primary_peak(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    working = utc_jja(frame)
    schedule = primary_peak_schedule(year)
    selected = working.merge(
        schedule, on=["analysis_date", "time"], how="inner",
        validate="many_to_one",
    )
    if selected.analysis_date.nunique() != 92 or \
            len(selected) != 92 * EXPECTED_SITES:
        raise ValueError(f"{year}: incomplete primary-grid peak fields")
    dense_means = (
        selected.groupby("analysis_date", observed=True).wbt.mean()
        .rename("regional_mean_wbt")
    )
    selected = selected.join(dense_means, on="analysis_date")
    return finalize(
        selected, year, "primary_grid_peak", selected.primary_regional_mean_wbt
    )


def process_year(path: Path, overwrite: bool) -> tuple[Path, dict]:
    import xarray as xr

    year = int(path.name.split("_")[2])
    output = OUTPUT_DIR / f"era5_land_{year}_jja_dense_daily_fields.csv.gz"
    with xr.open_dataset(path, engine="h5netcdf") as dataset:
        frame = dataset.to_dataframe().reset_index()
    frame["time"] = pd.to_datetime(frame.time, utc=True).dt.tz_localize(None)
    audit = validate_panel(frame, year)
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}")
        existing = pd.read_csv(
            output,
            usecols=["analysis_definition", "analysis_date", "time", "wbt"],
            parse_dates=["analysis_date", "time"],
        )
        peak_times = (
            existing[["analysis_definition", "analysis_date", "time"]]
            .drop_duplicates()
            .pivot(index="analysis_date", columns="analysis_definition",
                   values="time")
        )
        peak_shift_hours = (
            peak_times["dense_grid_peak"] - peak_times["primary_grid_peak"]
        ).dt.total_seconds() / 3600
        audit.update({
            "wbt_min_c": float(existing.wbt.min()),
            "wbt_max_c": float(existing.wbt.max()),
            "definitions": 2,
            "same_peak_days": int((peak_shift_hours == 0).sum()),
            "mean_absolute_peak_shift_hours": float(
                peak_shift_hours.abs().mean()
            ),
            "max_absolute_peak_shift_hours": float(
                peak_shift_hours.abs().max()
            ),
        })
        return output, audit

    frame["d2m"] = np.minimum(frame.d2m, frame.t2m)
    frame["wbt"] = PRIMARY.bolton_wbt_c(
        frame.t2m.to_numpy(), frame.d2m.to_numpy(), frame.sp.to_numpy()
    )
    if not np.isfinite(frame.wbt).all() or not frame.wbt.between(-60, 50).all():
        raise ValueError(f"{year}: invalid dense-grid WBT")

    fields = pd.concat([
        select_primary_peak(frame, year),
        select_dense_peak(frame, year),
    ], ignore_index=True)
    peak_times = (
        fields[["analysis_definition", "analysis_date", "time"]]
        .drop_duplicates()
        .pivot(index="analysis_date", columns="analysis_definition", values="time")
    )
    peak_shift_hours = (
        peak_times["dense_grid_peak"] - peak_times["primary_grid_peak"]
    ).dt.total_seconds() / 3600
    keep = [
        "analysis_role", "year", "month", "record_id", "analysis_date",
        "analysis_definition", "wbt_method", "time", "site_id",
        "requested_lon", "requested_lat", "label_mean_wbt",
        "regional_mean_wbt", "wbt", "t2m", "d2m", "sp",
    ]
    if "primary_regional_mean_wbt" in fields.columns:
        keep.append("primary_regional_mean_wbt")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fields[keep].to_csv(output, index=False, compression="gzip")
    audit["wbt_min_c"] = float(fields.wbt.min())
    audit["wbt_max_c"] = float(fields.wbt.max())
    audit["definitions"] = 2
    audit["same_peak_days"] = int((peak_shift_hours == 0).sum())
    audit["mean_absolute_peak_shift_hours"] = float(peak_shift_hours.abs().mean())
    audit["max_absolute_peak_shift_hours"] = float(peak_shift_hours.abs().max())
    print(f"[DONE] {output.name}: 2 definitions x 92 days x {EXPECTED_SITES} sites")
    return output, audit


def main() -> None:
    args = parse_args()
    paths = sorted(INPUT_DIR.glob("era5_land_*_jja_465sites.nc"))
    if args.years:
        requested = set(args.years)
        paths = [p for p in paths if int(p.name.split("_")[2]) in requested]
    if not paths:
        raise SystemExit(f"No dense yearly panels found in {INPUT_DIR}")
    audits = [process_year(path, args.overwrite)[1] for path in paths]
    AUDIT_FILE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audits).sort_values("year").to_csv(AUDIT_FILE, index=False)
    print(f"Audit: {AUDIT_FILE}")


if __name__ == "__main__":
    main()
