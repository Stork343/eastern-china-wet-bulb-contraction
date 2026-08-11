#!/usr/bin/env python3
"""Build frozen UTC daily fields for the 1950--1990 temporal extension."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
PROTOCOL_FILE = PROJECT_DIR / "EXTENSION_ANALYSIS_PROTOCOL.md"
SITE_FILE = PROJECT_DIR / "data" / "grid" / "eastern_china_121_sites.csv"
INPUT_DIR = PROJECT_DIR / "data" / "era5_historical_extension" / "hourly_points"
OUTPUT_DIR = PROJECT_DIR / "data" / "era5_historical_extension" / "daily_fields"
AUDIT_DIR = PROJECT_DIR / "output_historical_extension"
AUDIT_FILE = AUDIT_DIR / "historical_data_audit.csv"
MANIFEST_FILE = AUDIT_DIR / "historical_field_build_manifest.json"
YEARS = tuple(range(1950, 1991))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--years", type=int, nargs="+")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def validate_panel(frame: pd.DataFrame, year: int) -> dict[str, float | int | str]:
    required = {"time", "site_id", "requested_lon", "requested_lat",
                "t2m", "d2m", "sp"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{year}: missing columns {sorted(missing)}")
    if frame.time.nunique() != 2232 or frame.site_id.nunique() != 121:
        raise ValueError(f"{year}: incomplete buffered time-site panel")
    counts = frame.groupby("time", observed=True).site_id.nunique()
    if not (counts == 121).all():
        raise ValueError(f"{year}: incomplete hourly spatial field")
    for variable in ("t2m", "d2m", "sp"):
        if not np.isfinite(frame[variable]).all():
            raise ValueError(f"{year}: nonfinite values in {variable}")
    if not frame.t2m.between(220, 340).all():
        raise ValueError(f"{year}: implausible 2-m temperature")
    if not frame.sp.between(45000, 110000).all():
        raise ValueError(f"{year}: implausible surface pressure")
    return {
        "year": year,
        "analysis_role": "post-analysis historical extension",
        "forcing_segment": "1950-1978" if year <= 1978 else "1979-1990",
        "hours_with_buffer": frame.time.nunique(),
        "sites": frame.site_id.nunique(),
        "dewpoint_exceedances_clipped": int((frame.d2m > frame.t2m).sum()),
        "t2m_min_c": float(frame.t2m.min() - 273.15),
        "t2m_max_c": float(frame.t2m.max() - 273.15),
        "sp_min_pa": float(frame.sp.min()),
        "sp_max_pa": float(frame.sp.max()),
    }


def select_utc_peak_fields(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    working = frame.copy()
    working["analysis_date"] = working.time.dt.floor("D")
    working = working.loc[working.analysis_date.dt.month.isin((6, 7, 8))].copy()
    regional = (working.groupby("time", observed=True).wbt.mean()
                .rename("regional_mean_wbt").reset_index())
    regional["analysis_date"] = regional.time.dt.floor("D")
    peaks = (regional.sort_values(
        ["analysis_date", "regional_mean_wbt", "time"],
        ascending=[True, False, True]
    ).drop_duplicates("analysis_date"))
    selected = working.merge(
        peaks[["analysis_date", "time", "regional_mean_wbt"]],
        on=["analysis_date", "time"], how="inner", validate="many_to_one"
    )
    if selected.analysis_date.nunique() != 92 or len(selected) != 92 * 121:
        raise ValueError(f"{year}: incomplete UTC daily peak fields")
    selected["year"] = year
    selected["month"] = selected.analysis_date.dt.month
    selected["record_id"] = selected.analysis_date.dt.strftime("%Y%m")
    selected["day_definition"] = "utc"
    selected["wbt_method"] = "pressure_aware_Bolton"
    selected["analysis_role"] = "post-analysis historical extension"
    selected["forcing_segment"] = "1950-1978" if year <= 1978 else "1979-1990"
    return selected


def process_year(path: Path, overwrite: bool) -> tuple[Path, dict]:
    import xarray as xr

    year = int(path.name.split("_")[2])
    output = OUTPUT_DIR / f"era5_land_{year}_jja_daily_fields.csv.gz"
    with xr.open_dataset(path, engine="h5netcdf") as dataset:
        frame = dataset.to_dataframe().reset_index()
    frame["time"] = pd.to_datetime(frame.time, utc=True).dt.tz_localize(None)
    audit = validate_panel(frame, year)
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}")
        return output, audit

    frame["d2m"] = np.minimum(frame.d2m, frame.t2m)
    frame["wbt"] = bolton_wbt_c(
        frame.t2m.to_numpy(), frame.d2m.to_numpy(), frame.sp.to_numpy()
    )
    if not np.isfinite(frame.wbt).all() or not frame.wbt.between(-60, 50).all():
        raise ValueError(f"{year}: invalid wet-bulb temperatures")
    fields = select_utc_peak_fields(frame, year)
    keep = [
        "analysis_role", "forcing_segment", "year", "month", "record_id",
        "analysis_date", "day_definition", "wbt_method", "time", "site_id",
        "requested_lon", "requested_lat", "regional_mean_wbt", "wbt",
        "t2m", "d2m", "sp",
    ]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fields[keep].to_csv(output, index=False, compression="gzip")
    audit["wbt_min_c"] = float(fields.wbt.min())
    audit["wbt_max_c"] = float(fields.wbt.max())
    audit["input_sha256"] = sha256(path)
    audit["output_sha256"] = sha256(output)
    print(f"[DONE] {output.name}: 92 UTC peak fields x 121 sites")
    return output, audit


def main() -> None:
    args = parse_args()
    paths = sorted(INPUT_DIR.glob("era5_land_*_jja_121sites.nc"))
    requested_years = set(args.years or YEARS)
    unknown = sorted(requested_years - set(YEARS))
    if unknown:
        raise SystemExit(f"Years outside frozen 1950-1990 period: {unknown}")
    paths = [p for p in paths if int(p.name.split("_")[2]) in requested_years]
    found = {int(p.name.split("_")[2]) for p in paths}
    missing = sorted(requested_years - found)
    if missing:
        raise SystemExit(f"Missing yearly point panels: {missing}")
    results = [process_year(path, args.overwrite) for path in paths]
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([audit for _, audit in results]).sort_values("year").to_csv(
        AUDIT_FILE, index=False
    )
    outputs = [path for path, _ in results]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis historical extension",
        "protocol_sha256": sha256(PROTOCOL_FILE),
        "site_manifest_sha256": sha256(SITE_FILE),
        "years": sorted(requested_years),
        "field_definition": "UTC hour maximizing the 121-site regional mean",
        "wet_bulb_method": "pressure-aware Bolton implementation from script 22",
        "outputs": [
            {"file": str(path.relative_to(PROJECT_DIR)), "sha256": sha256(path)}
            for path in outputs
        ],
        "audit_file": str(AUDIT_FILE.relative_to(PROJECT_DIR)),
    }
    MANIFEST_FILE.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Audit: {AUDIT_FILE}")
    print(f"Manifest: {MANIFEST_FILE}")


if __name__ == "__main__":
    main()
