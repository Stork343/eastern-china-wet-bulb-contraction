#!/usr/bin/env python3
"""Download the frozen non-development NOAA ISD station-year panel."""

from __future__ import annotations

import hashlib
import json
import math
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent.parent
PROTOCOL = PROJECT / "EXTENSION_ANALYSIS_PROTOCOL.md"
SITE_FILE = PROJECT / "data" / "grid" / "eastern_china_121_sites.csv"
BASE = PROJECT / "data" / "noaa_isd_extension"
RAW = BASE / "raw_global_hourly"
HISTORY_FILE = BASE / "isd-history_20260809.csv"
MANIFEST_FILE = BASE / "noaa_extension_station_manifest.csv"
QUALIFICATION_FILE = BASE / "noaa_extension_station_year_qualification.csv"
PROVENANCE_FILE = BASE / "noaa_extension_download_audit.json"
HISTORY_URL = "https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv"
DATA_URL = "https://www.ncei.noaa.gov/data/global-hourly/access/{year}/{station}.csv"
YEARS = (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023, 2025)
REQUESTED_END = 20250831
STATIONS = 30
MIN_COMPLETE_HOURS = 400
ACCEPTED_QC = {"0", "1", "4", "5", "9", "A", "C", "I", "M", "P", "R", "U"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path, attempts: int = 8) -> None:
    if path.exists() and path.stat().st_size > 100:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    import requests

    for attempt in range(1, attempts + 1):
        try:
            with requests.get(
                url, headers={"User-Agent": "JRSSC-reproducibility/1.0"},
                timeout=(30, 120), stream=True
            ) as response:
                response.raise_for_status()
                with partial.open("wb") as output:
                    for block in response.iter_content(1024 * 1024):
                        if block:
                            output.write(block)
            os.replace(partial, path)
            return
        except (requests.RequestException, TimeoutError, ConnectionError):
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(min(2**attempt, 30))


def distance_km(lon1, lat1, lon2, lat2):
    mean_lat = np.deg2rad((lat1 + lat2) / 2)
    dx = (lon1 - lon2) * 111.32 * np.cos(mean_lat)
    dy = (lat1 - lat2) * 110.57
    return np.sqrt(dx * dx + dy * dy)


def candidate_panel(history: pd.DataFrame, sites: pd.DataFrame, end_rule: int) -> pd.DataFrame:
    first = min(YEARS) * 10000 + 601
    frame = history.loc[
        history.LAT.between(20, 42) & history.LON.between(105, 125) &
        history.BEGIN.le(first) & history.END.ge(end_rule)
    ].copy()
    if frame.empty:
        return frame
    frame["nearest_analysis_site_km"] = [
        float(np.min(distance_km(row.LON, row.LAT, sites.lon, sites.lat)))
        for row in frame.itertuples(index=False)
    ]
    frame = frame.loc[frame.nearest_analysis_site_km.le(150)].copy()
    frame["USAF"] = frame.USAF.astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(6)
    frame["WBAN"] = frame.WBAN.astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)
    frame["station"] = frame.USAF + frame.WBAN
    return frame.drop_duplicates("station").sort_values("station").reset_index(drop=True)


def maximin(candidates: pd.DataFrame, count: int) -> pd.DataFrame:
    if len(candidates) < count:
        raise ValueError(f"Only {len(candidates)} candidates satisfy the administrative rules")
    center = distance_km(candidates.LON, candidates.LAT, 115, 31)
    selected = [int(np.argmin(center))]
    while len(selected) < count:
        remaining = [i for i in candidates.index if i not in selected]
        score = []
        for i in remaining:
            score.append(min(
                float(distance_km(candidates.at[i, "LON"], candidates.at[i, "LAT"],
                                  candidates.at[j, "LON"], candidates.at[j, "LAT"]))
                for j in selected
            ))
        best_score = max(score)
        # Ties are broken by the lexicographically sorted station order.
        best = min(i for i, value in zip(remaining, score)
                   if math.isclose(value, best_score, rel_tol=0, abs_tol=1e-12))
        selected.append(best)
    result = candidates.loc[selected].copy().reset_index(drop=True)
    result.insert(0, "selection_order", np.arange(1, len(result) + 1))
    return result


def parse_isd(series: pd.Series, missing: int) -> tuple[pd.Series, pd.Series]:
    parts = series.fillna("").astype(str).str.split(",", n=2, expand=True)
    values = pd.to_numeric(parts[0], errors="coerce")
    values = values.mask(values.abs().eq(missing)) / 10
    quality = parts[1].astype(str) if parts.shape[1] > 1 else pd.Series("", index=series.index)
    return values, quality


def station_pressure(slp_hpa: pd.Series, elevation_m: pd.Series) -> pd.Series:
    factor = np.maximum(1 - 2.25577e-5 * elevation_m, 0.1) ** 5.2559
    return slp_hpa * 100 * factor


def process_source(path: Path, station: str, year: int) -> tuple[pd.DataFrame, dict]:
    columns = ["STATION", "DATE", "LATITUDE", "LONGITUDE", "ELEVATION", "NAME",
               "REPORT_TYPE", "TMP", "DEW", "SLP"]
    frame = pd.read_csv(path, usecols=lambda name: name in columns, dtype={"STATION": str},
                        low_memory=False)
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{station} {year} lacks columns {sorted(missing)}")
    frame["time_utc"] = pd.to_datetime(frame.DATE, errors="coerce", utc=True)
    frame = frame.loc[
        frame.time_utc.dt.month.isin((6, 7, 8)) &
        frame.time_utc.dt.minute.eq(0) & frame.time_utc.dt.second.eq(0)
    ].copy()
    frame["temperature_c"], frame["temperature_qc"] = parse_isd(frame.TMP, 9999)
    frame["dewpoint_c"], frame["dewpoint_qc"] = parse_isd(frame.DEW, 9999)
    frame["slp_hpa"], frame["slp_qc"] = parse_isd(frame.SLP, 99999)
    frame["station_pressure_pa"] = station_pressure(frame.slp_hpa, frame.ELEVATION)
    frame["qc_accepted"] = (
        frame.temperature_qc.isin(ACCEPTED_QC) &
        frame.dewpoint_qc.isin(ACCEPTED_QC) & frame.slp_qc.isin(ACCEPTED_QC)
    )
    frame["complete"] = (
        frame[["temperature_c", "dewpoint_c", "station_pressure_pa"]].notna().all(axis=1) &
        frame.qc_accepted & frame.station_pressure_pa.between(45000, 110000) &
        frame.temperature_c.between(-60, 55) & frame.dewpoint_c.between(-80, 40)
    )
    frame = frame.sort_values(["time_utc", "REPORT_TYPE"]).drop_duplicates("time_utc")
    complete = frame.loc[frame.complete].copy()
    complete["dewpoint_c"] = np.minimum(complete.dewpoint_c, complete.temperature_c)
    complete["station"] = station
    complete["year"] = year
    keep = ["station", "year", "time_utc", "LATITUDE", "LONGITUDE", "ELEVATION", "NAME",
            "REPORT_TYPE", "temperature_c", "temperature_qc", "dewpoint_c", "dewpoint_qc",
            "slp_hpa", "slp_qc", "station_pressure_pa"]
    summary = {
        "station": station,
        "year": year,
        "jja_exact_hour_rows": int(len(frame)),
        "complete_exact_hours": int(len(complete)),
        "qualified": bool(len(complete) >= MIN_COMPLETE_HOURS),
        "source_file": str(path.relative_to(PROJECT)),
        "source_sha256": sha256(path),
    }
    return complete[keep], summary


def main() -> None:
    BASE.mkdir(parents=True, exist_ok=True)
    download(HISTORY_URL, HISTORY_FILE)
    history = pd.read_csv(HISTORY_FILE, dtype={"USAF": str, "WBAN": str})
    history["BEGIN"] = pd.to_numeric(history.BEGIN, errors="coerce")
    history["END"] = pd.to_numeric(history.END, errors="coerce")
    sites = pd.read_csv(SITE_FILE)
    snapshot_global_max_end = int(history.END.max())
    history_in_scope = history.loc[
        history.LAT.between(20, 42) & history.LON.between(105, 125) &
        history.BEGIN.le(min(YEARS) * 10000 + 601)
    ]
    snapshot_max_end = int(history_in_scope.END.max())
    original = candidate_panel(history, sites, REQUESTED_END)
    adjusted = candidate_panel(history, sites, snapshot_max_end)
    if len(original) != 0:
        raise RuntimeError(
            "The frozen infeasibility audit changed: END>=20250831 now has candidates"
        )
    stations = maximin(adjusted, STATIONS)
    stations["selection_rule"] = "deterministic maximin; no outcome values used"
    stations["administrative_end_rule"] = snapshot_max_end
    stations.to_csv(MANIFEST_FILE, index=False)

    summaries: list[dict] = []
    yearly_outputs = []
    errors = []
    for year in YEARS:
        pieces = []

        def one_station(row):
            station = str(row.station)
            url = DATA_URL.format(year=year, station=station)
            raw_path = RAW / str(year) / f"{station}.csv"
            try:
                download(url, raw_path)
                piece, summary = process_source(raw_path, station, year)
                summary["source_url"] = url
                return piece, summary, None
            except Exception as error:  # retained in the machine-readable audit
                item = {"station": station, "year": year, "url": url,
                        "error": f"{type(error).__name__}: {error}"}
                summary = {"station": station, "year": year,
                           "jja_exact_hour_rows": 0, "complete_exact_hours": 0,
                           "qualified": False, "source_url": url,
                           "error": item["error"]}
                return None, summary, item

        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = [executor.submit(one_station, row)
                       for row in stations.itertuples(index=False)]
            for future in as_completed(futures):
                piece, summary, error = future.result()
                summaries.append(summary)
                if piece is not None:
                    pieces.append(piece)
                if error is not None:
                    errors.append(error)
        combined = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
        output = BASE / f"noaa_isd_extension_{year}_jja_exact_hours.csv.gz"
        combined.to_csv(output, index=False, compression="gzip")
        yearly_outputs.append(output)
        print(f"{year}: {combined.station.nunique() if len(combined) else 0} stations, "
              f"{len(combined)} complete exact hours", flush=True)

    qualification = pd.DataFrame(summaries).sort_values(["year", "station"])
    qualification.to_csv(QUALIFICATION_FILE, index=False)
    if qualification.groupby("year").qualified.sum().min() < 1:
        raise RuntimeError("At least one evaluation year has no qualified station-year")
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis non-development station extension",
        "protocol_file": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "history_url": HISTORY_URL,
        "history_file": str(HISTORY_FILE.relative_to(PROJECT)),
        "history_sha256": sha256(HISTORY_FILE),
        "history_snapshot_global_max_end": snapshot_global_max_end,
        "history_snapshot_in_scope_max_end": snapshot_max_end,
        "frozen_requested_end": REQUESTED_END,
        "original_candidate_count": len(original),
        "operational_candidate_count": len(adjusted),
        "operationalization": (
            "The official history snapshot ended before 2025-08-31; within the "
            "fixed study rectangle its latest END was 2025-08-24, so the frozen "
            "END>=2025-08-31 rule yielded zero candidates. The administrative "
            "end filter was truncated to the in-scope snapshot maximum; years, geography, "
            "maximin count and outcome-blind selection were unchanged."
        ),
        "years": list(YEARS),
        "selected_stations": len(stations),
        "minimum_complete_exact_hours": MIN_COMPLETE_HOURS,
        "selection_uses_effects": False,
        "errors": errors,
        "outputs": [
            {"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
            for path in (MANIFEST_FILE, QUALIFICATION_FILE, *yearly_outputs)
        ],
    }
    PROVENANCE_FILE.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(f"Selected {len(stations)} stations from {len(adjusted)} operational candidates")
    print(qualification.groupby("year").agg(
        retrieved=("station", "count"), qualified=("qualified", "sum"),
        complete_hours=("complete_exact_hours", "sum")
    ).to_string())


if __name__ == "__main__":
    main()
