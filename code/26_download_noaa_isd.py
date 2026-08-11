#!/usr/bin/env python3
"""Acquire an independent NOAA ISD hourly validation sample.

Stations must lie inside the study domain and span both discovery summers.
A deterministic maximin rule selects a geographically dispersed subset. Only
JJA temperature, dew point, sea-level pressure, coordinates, and QC flags are
retained from the public yearly CSV files.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_DIR / "data" / "noaa_isd_validation"
SITE_FILE = PROJECT_DIR / "data" / "grid" / "eastern_china_121_sites.csv"
HISTORY_URL = "https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv"
DATA_URL = "https://www.ncei.noaa.gov/data/global-hourly/access/{year}/{station}.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--years", type=int, nargs="+", default=[2015, 2022])
    parser.add_argument("--stations", type=int, default=30)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def distance_km(lon1, lat1, lon2, lat2):
    mean_lat = np.deg2rad((lat1 + lat2) / 2)
    dx = (lon1 - lon2) * 111.32 * np.cos(mean_lat)
    dy = (lat1 - lat2) * 110.57
    return np.sqrt(dx * dx + dy * dy)


def select_stations(history: pd.DataFrame, sites: pd.DataFrame,
                    count: int, years: list[int]) -> pd.DataFrame:
    first = min(years) * 10000 + 601
    last = max(years) * 10000 + 831
    candidates = history.loc[
        history["LAT"].between(20, 42) & history["LON"].between(105, 125) &
        (history["BEGIN"] <= first) & (history["END"] >= last)
    ].copy()
    candidates["nearest_analysis_site_km"] = [
        float(np.min(distance_km(row.LON, row.LAT, sites["lon"], sites["lat"])))
        for row in candidates.itertuples(index=False)
    ]
    candidates = candidates.loc[candidates["nearest_analysis_site_km"] <= 150].copy()
    candidates["station"] = (
        candidates["USAF"].astype(str).str.zfill(6) +
        candidates["WBAN"].astype(str).str.zfill(5)
    )
    candidates = candidates.drop_duplicates("station").reset_index(drop=True)
    if len(candidates) < count:
        raise ValueError(f"Only {len(candidates)} stations span the requested years")

    center_distance = distance_km(candidates["LON"], candidates["LAT"], 115, 31)
    selected = [int(np.argmin(center_distance))]
    while len(selected) < count:
        remaining = [i for i in candidates.index if i not in selected]
        minimum_distance = []
        for i in remaining:
            distances = [distance_km(candidates.at[i, "LON"], candidates.at[i, "LAT"],
                                     candidates.at[j, "LON"], candidates.at[j, "LAT"])
                         for j in selected]
            minimum_distance.append(min(distances))
        selected.append(remaining[int(np.argmax(minimum_distance))])
    return candidates.loc[selected].reset_index(drop=True)


def parse_isd_value(series: pd.Series, missing: int = 9999) -> tuple[pd.Series, pd.Series]:
    parts = series.fillna("").str.split(",", n=2, expand=True)
    value = pd.to_numeric(parts[0], errors="coerce")
    quality = parts[1] if parts.shape[1] > 1 else pd.Series(index=series.index, dtype=str)
    value = value.where(value.abs() != missing)
    return value / 10, quality


def standard_pressure_from_sea_level(slp_hpa: pd.Series,
                                     elevation_m: pd.Series) -> pd.Series:
    factor = np.maximum(1 - 2.25577e-5 * elevation_m, 0.1) ** 5.2559
    return slp_hpa * 100 * factor


def download_year(stations: pd.DataFrame, year: int, overwrite: bool) -> pd.DataFrame:
    output = OUTPUT_DIR / f"noaa_isd_{year}_jja.csv.gz"
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}")
        return pd.read_csv(output)

    pieces = []
    columns = ["STATION", "DATE", "LATITUDE", "LONGITUDE", "ELEVATION",
               "NAME", "REPORT_TYPE", "QUALITY_CONTROL", "TMP", "DEW", "SLP"]
    for row in stations.itertuples(index=False):
        url = DATA_URL.format(year=year, station=row.station)
        try:
            frame = pd.read_csv(url, usecols=columns, low_memory=False)
        except Exception as error:
            print(f"[WARN] {row.station} {year}: {error}")
            continue
        frame["DATE"] = pd.to_datetime(frame["DATE"], errors="coerce", utc=True)
        frame = frame.loc[frame["DATE"].dt.month.isin((6, 7, 8))].copy()
        if frame.empty:
            continue
        frame["temperature_c"], frame["temperature_qc"] = parse_isd_value(frame["TMP"])
        frame["dewpoint_c"], frame["dewpoint_qc"] = parse_isd_value(frame["DEW"])
        frame["slp_hpa"], frame["slp_qc"] = parse_isd_value(frame["SLP"], missing=99999)
        frame["station_pressure_pa"] = standard_pressure_from_sea_level(
            frame["slp_hpa"], frame["ELEVATION"]
        )
        frame["source_url"] = url
        pieces.append(frame[[
            "STATION", "DATE", "LATITUDE", "LONGITUDE", "ELEVATION", "NAME",
            "REPORT_TYPE", "QUALITY_CONTROL", "temperature_c", "temperature_qc",
            "dewpoint_c", "dewpoint_qc", "slp_hpa", "slp_qc",
            "station_pressure_pa", "source_url",
        ]])
        print(f"[DONE] {row.station} {year}: {len(frame)} JJA observations")
    if not pieces:
        raise RuntimeError(f"No NOAA observations retrieved for {year}")
    combined = pd.concat(pieces, ignore_index=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False, compression="gzip")
    return combined


def main() -> None:
    args = parse_args()
    history = pd.read_csv(HISTORY_URL)
    history["BEGIN"] = pd.to_numeric(history["BEGIN"], errors="coerce")
    history["END"] = pd.to_numeric(history["END"], errors="coerce")
    sites = pd.read_csv(SITE_FILE)
    stations = select_stations(history, sites, args.stations, args.years)
    print(stations[["station", "STATION NAME", "CTRY", "LAT", "LON", "ELEV(M)",
                    "nearest_analysis_site_km"]]
          .to_string(index=False))
    if args.dry_run:
        print("Dry run complete; no yearly station files were downloaded.")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stations.to_csv(OUTPUT_DIR / "noaa_isd_station_manifest.csv", index=False)
    summaries = []
    for year in args.years:
        frame = download_year(stations, year, args.overwrite)
        summaries.append({
            "year": year,
            "stations_retrieved": int(frame["STATION"].nunique()),
            "observations": len(frame),
            "complete_temperature_dewpoint_pressure": int(frame[
                ["temperature_c", "dewpoint_c", "station_pressure_pa"]
            ].notna().all(axis=1).sum()),
        })
    (OUTPUT_DIR / "noaa_isd_provenance.json").write_text(json.dumps({
        "history_url": HISTORY_URL,
        "yearly_url_template": DATA_URL,
        "selection": "deterministic maximin within the study domain",
        "years": args.years,
        "station_count_requested": args.stations,
        "summaries": summaries,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
