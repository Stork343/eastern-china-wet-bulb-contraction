#!/usr/bin/env python3
"""Acquire ERA5-Land point series at the frozen NOAA extension stations."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent.parent
PROTOCOL = PROJECT / "EXTENSION_ANALYSIS_PROTOCOL.md"
MANIFEST = PROJECT / "data" / "noaa_isd_extension" / "noaa_extension_station_manifest.csv"
BASE = PROJECT / "data" / "noaa_isd_extension" / "era5_land_points"
RAW = BASE / "raw_archives"
TRIMMED = BASE / "trimmed_points"
PROVENANCE = BASE / "era5_station_point_provenance.json"
DATASET = "reanalysis-era5-land-timeseries"
DATE_RANGE = "1992-05-31/2025-09-01"
YEARS = (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023, 2025)
VARIABLES = ("2m_dewpoint_temperature", "2m_temperature", "surface_pressure")
EXPECTED = {"d2m", "t2m", "sp"}


class NoFiniteLandCell(RuntimeError):
    """The nearest ERA5-Land grid cell is outside the product's land mask."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def requested_index() -> pd.DatetimeIndex:
    pieces = [pd.date_range(f"{year}-05-31 16:00:00", f"{year}-09-01 15:00:00", freq="h")
              for year in YEARS]
    return pieces[0].append(pieces[1:])


def valid_zip(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 1000 and zipfile.is_zipfile(path)


def client():
    import cdsapi
    try:
        return cdsapi.Client(quiet=True, progress=False, maximum_tries=3,
                             retry_after=30, sleep_max=30)
    except TypeError as error:
        if "unexpected keyword argument" not in str(error):
            raise
        return cdsapi.Client(quiet=True, progress=False, retry_max=3,
                             sleep_max=30)


def acquire(row, attempts: int) -> Path:
    station = str(row.station)
    output = RAW / f"station_{station}.zip"
    if valid_zip(output):
        return output
    partial = output.with_suffix(".zip.part")
    request = {
        "variable": list(VARIABLES),
        "location": {"latitude": float(row.LAT), "longitude": float(row.LON)},
        "date": [DATE_RANGE], "data_format": "netcdf",
    }
    for attempt in range(1, attempts + 1):
        partial.unlink(missing_ok=True)
        try:
            client().retrieve(DATASET, request, str(partial))
            if not valid_zip(partial):
                raise RuntimeError("CDS response is not a valid ZIP archive")
            os.replace(partial, output)
            return output
        except Exception as error:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            print(
                f"[RETRY] station {station}, attempt {attempt}/{attempts}: "
                f"{type(error).__name__}: {error}",
                flush=True,
            )
            time.sleep(min(15 * attempt, 60))
    raise AssertionError("unreachable")


def trim(row, archive: Path) -> Path:
    import xarray as xr

    station = str(row.station)
    output = TRIMMED / f"station_{station}_evaluation_jja.nc"
    if output.exists() and output.stat().st_size > 1000:
        return output
    with tempfile.TemporaryDirectory(prefix=f"era5-station-{station}-") as directory:
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(directory)
        files = sorted(Path(directory).rglob("*.nc"))
        datasets = []
        for path in files:
            with xr.open_dataset(path, engine="h5netcdf") as source:
                datasets.append(source.load())
        panel = xr.merge(datasets, compat="override", join="exact")
    if "valid_time" in panel.coords:
        panel = panel.rename({"valid_time": "time"})
    if set(panel.data_vars) != EXPECTED:
        raise RuntimeError(f"{station}: variables {set(panel.data_vars)}")
    index = requested_index()
    positions = pd.DatetimeIndex(panel.time.values).get_indexer(index)
    if (positions < 0).any():
        raise RuntimeError(f"{station}: missing {(positions < 0).sum()} requested hours")
    panel = panel.isel(time=positions)[sorted(EXPECTED)].assign_coords(time=index)
    if panel.sizes["time"] != 2232 * len(YEARS):
        raise RuntimeError(f"{station}: incorrect trimmed time count")
    if all(not np.isfinite(panel[variable]).any() for variable in EXPECTED):
        raise NoFiniteLandCell(
            f"nearest grid cell ({float(panel.latitude.values):.1f}, "
            f"{float(panel.longitude.values):.1f}) has no finite ERA5-Land values"
        )
    for variable in EXPECTED:
        if not np.isfinite(panel[variable]).all():
            raise RuntimeError(f"{station}: nonfinite {variable}")
    panel.attrs.update({
        "station": station, "requested_lat": float(row.LAT),
        "requested_lon": float(row.LON), "source_dataset": DATASET,
        "source_archive": archive.name, "protocol_sha256": sha256(PROTOCOL),
    })
    encoding = {name: {"zlib": True, "complevel": 4, "shuffle": True}
                for name in panel.data_vars}
    partial = output.with_suffix(".nc.part")
    panel.to_netcdf(partial, engine="h5netcdf", encoding=encoding)
    os.replace(partial, output)
    return output


def one(row, attempts):
    archive = acquire(row, attempts)
    try:
        output = trim(row, archive)
    except NoFiniteLandCell as error:
        print(f"[UNAVAILABLE] {row.station}: {error}", flush=True)
        return archive, None, str(error)
    print(f"[DONE] {row.station}", flush=True)
    return archive, output, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--attempts", type=int, default=12)
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    TRIMMED.mkdir(parents=True, exist_ok=True)
    stations = pd.read_csv(MANIFEST, dtype={"station": str}).sort_values("selection_order")
    if len(stations) != 30:
        raise RuntimeError("Expected the frozen 30-station manifest")
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(one, row, args.attempts): str(row.station)
                   for row in stations.itertuples(index=False)}
        for future in as_completed(futures):
            results.append((futures[future], *future.result()))
    results.sort()
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "non-development station measurement extension",
        "dataset": DATASET, "date_range_requested": DATE_RANGE,
        "years_retained": list(YEARS), "variables": list(VARIABLES),
        "protocol_sha256": sha256(PROTOCOL), "station_manifest_sha256": sha256(MANIFEST),
        "selected_stations": len(results),
        "available_station_series": sum(output is not None
                                         for _, _, output, _ in results),
        "stations": [
            {
                "station": station,
                "archive": str(archive.relative_to(PROJECT)),
                "archive_sha256": sha256(archive),
                "status": "available" if output is not None else
                          "unavailable_nearest_land_mask_cell",
                "trimmed": str(output.relative_to(PROJECT))
                           if output is not None else None,
                "trimmed_sha256": sha256(output) if output is not None else None,
                "reason": reason,
            }
            for station, archive, output, reason in results
        ],
    }
    PROVENANCE.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(
        f"Completed {audit['available_station_series']} finite ERA5-Land "
        f"station point series from {len(results)} frozen stations"
    )


if __name__ == "__main__":
    main()
