#!/usr/bin/env python3
"""Acquire the frozen 1950--1990 ERA5-Land point-panel extension.

This post-analysis extension is isolated from the 1991--2025 confirmatory
cache.  The CDS time-series product is queried once per frozen site for the
three variables needed by the pressure-aware Bolton wet-bulb calculation.
Raw archives, trimmed JJA buffers, yearly panels, hashes, and request metadata
are retained so interrupted runs resume without changing the design.
"""

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


PROJECT_DIR = Path(__file__).resolve().parent.parent
PROTOCOL_FILE = PROJECT_DIR / "EXTENSION_ANALYSIS_PROTOCOL.md"
SITE_FILE = PROJECT_DIR / "data" / "grid" / "eastern_china_121_sites.csv"
BASE_DIR = PROJECT_DIR / "data" / "era5_historical_extension"
RAW_DIR = BASE_DIR / "cds_point_archives"
SLIM_DIR = BASE_DIR / "trimmed_points"
OUTPUT_DIR = BASE_DIR / "hourly_points"
PROVENANCE_FILE = BASE_DIR / "historical_cds_point_provenance.json"

DATASET = "reanalysis-era5-land-timeseries"
ANALYSIS_YEARS = tuple(range(1950, 1991))
DATE_RANGE = "1950-05-31/1990-09-01"
VARIABLES = (
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_pressure",
)
EXPECTED_NAMES = {"d2m", "t2m", "sp"}
EXPECTED_HOURS_PER_YEAR = 2232


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-attempts", type=int, default=12)
    parser.add_argument("--sites", type=int, nargs="+")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--assemble-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid_zip(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 1000 and zipfile.is_zipfile(path)


def requested_time_index() -> pd.DatetimeIndex:
    pieces = [
        pd.date_range(
            f"{year}-05-31 16:00:00",
            f"{year}-09-01 15:00:00",
            freq="h",
        )
        for year in ANALYSIS_YEARS
    ]
    return pieces[0].append(pieces[1:])


def make_client():
    import cdsapi

    try:
        return cdsapi.Client(
            quiet=True, progress=False, maximum_tries=3,
            retry_after=30, sleep_max=30,
        )
    except TypeError as error:
        if "unexpected keyword argument" not in str(error):
            raise
        return cdsapi.Client(
            quiet=True, progress=False, retry_max=3, sleep_max=30
        )


def download_archive(site: pd.Series, overwrite: bool,
                     max_attempts: int) -> Path:
    site_id = int(site.site_id)
    output = RAW_DIR / f"site_{site_id:03d}.zip"
    if valid_zip(output) and not overwrite:
        return output
    partial = output.with_suffix(".zip.part")
    partial.unlink(missing_ok=True)
    request = {
        "variable": list(VARIABLES),
        "location": {
            "latitude": float(site.lat),
            "longitude": float(site.lon),
        },
        "date": [DATE_RANGE],
        "data_format": "netcdf",
    }
    for attempt in range(1, max_attempts + 1):
        try:
            make_client().retrieve(DATASET, request, str(partial))
            break
        except Exception as error:
            partial.unlink(missing_ok=True)
            if attempt == max_attempts:
                raise
            delay = min(15 * attempt, 60)
            print(
                f"[RETRY] site {site_id:03d}, attempt {attempt}/{max_attempts}: "
                f"{type(error).__name__}: {error}; waiting {delay}s",
                flush=True,
            )
            time.sleep(delay)
    if not valid_zip(partial):
        raise ValueError(f"Invalid CDS archive for site {site_id}")
    os.replace(partial, output)
    return output


def trim_archive(site: pd.Series, archive: Path, overwrite: bool) -> Path:
    import xarray as xr

    site_id = int(site.site_id)
    output = SLIM_DIR / f"site_{site_id:03d}_jja_buffers.nc"
    if output.exists() and not overwrite:
        return output

    with tempfile.TemporaryDirectory(prefix=f"era5-historical-{site_id:03d}-") as tmp:
        extraction = Path(tmp)
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(extraction)
        files = sorted(extraction.rglob("*.nc"))
        if not files:
            raise ValueError(f"No NetCDF members in {archive}")
        datasets = []
        for path in files:
            with xr.open_dataset(path, engine="h5netcdf") as source:
                datasets.append(source.load())
        panel = xr.merge(datasets, compat="override", join="exact")

    if "valid_time" in panel.coords:
        panel = panel.rename({"valid_time": "time"})
    missing = EXPECTED_NAMES - set(panel.data_vars)
    if missing:
        raise ValueError(f"Site {site_id} lacks variables {sorted(missing)}")
    actual_lat = float(panel.coords["latitude"].values)
    actual_lon = float(panel.coords["longitude"].values)
    if (abs(actual_lat - float(site.lat)) > 0.051 or
            abs(actual_lon - float(site.lon)) > 0.051):
        raise ValueError(f"Unexpected matched coordinate for site {site_id}")

    expected_time = requested_time_index()
    available = pd.DatetimeIndex(panel.time.values)
    positions = available.get_indexer(expected_time)
    if (positions < 0).any():
        raise ValueError(f"Site {site_id} lacks {int((positions < 0).sum())} hours")
    panel = panel.isel(time=positions)[sorted(EXPECTED_NAMES)]
    expected_count = EXPECTED_HOURS_PER_YEAR * len(ANALYSIS_YEARS)
    if panel.sizes.get("time") != expected_count:
        raise ValueError(f"Unexpected trimmed time count for site {site_id}")
    for name in EXPECTED_NAMES:
        if not np.isfinite(panel[name].values).all():
            raise ValueError(f"Nonfinite {name} values at site {site_id}")

    panel = panel.assign_coords(time=expected_time)
    panel.attrs.update({
        "site_id": site_id,
        "requested_lon": float(site.lon),
        "requested_lat": float(site.lat),
        "matched_lon": actual_lon,
        "matched_lat": actual_lat,
        "source_archive": archive.name,
        "source_dataset": DATASET,
        "analysis_period": "1950-1990 post-analysis temporal extension",
        "protocol_sha256": sha256(PROTOCOL_FILE),
    })
    encoding = {
        name: {"zlib": True, "complevel": 4, "shuffle": True}
        for name in panel.data_vars
    }
    partial = output.with_suffix(".nc.part")
    partial.unlink(missing_ok=True)
    panel.to_netcdf(partial, engine="h5netcdf", encoding=encoding)
    os.replace(partial, output)
    return output


def process_site(site: pd.Series, overwrite: bool, assemble_only: bool,
                 max_attempts: int) -> Path:
    site_id = int(site.site_id)
    archive = RAW_DIR / f"site_{site_id:03d}.zip"
    if assemble_only:
        if not valid_zip(archive):
            raise FileNotFoundError(f"Missing archive for site {site_id}")
    else:
        archive = download_archive(site, overwrite, max_attempts)
    output = trim_archive(site, archive, overwrite)
    print(f"[DONE] site {site_id:03d}", flush=True)
    return output


def assemble_year(site_files: list[Path], sites: pd.DataFrame, year: int,
                  overwrite: bool) -> Path:
    import xarray as xr

    output = OUTPUT_DIR / f"era5_land_{year}_jja_121sites.nc"
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}", flush=True)
        return output
    start = f"{year}-05-31 16:00:00"
    end = f"{year}-09-01 15:00:00"
    pieces = []
    for path, site in zip(site_files, sites.itertuples(index=False)):
        with xr.open_dataset(path, engine="h5netcdf") as source:
            piece = source.sel(time=slice(start, end)).load()
        pieces.append(piece.expand_dims(site=[int(site.site_id)]))
    panel = xr.concat(pieces, dim="site", coords="minimal", compat="override")
    panel = panel.assign_coords(
        site_id=("site", sites.site_id.to_numpy(dtype=np.int32)),
        requested_lon=("site", sites.lon.to_numpy()),
        requested_lat=("site", sites.lat.to_numpy()),
    )
    if (panel.sizes.get("time") != EXPECTED_HOURS_PER_YEAR or
            panel.sizes.get("site") != 121):
        raise ValueError(f"Unexpected dimensions for {year}: {dict(panel.sizes)}")
    panel.attrs.update({
        "spatial_manifest": SITE_FILE.name,
        "spatial_manifest_sha256": sha256(SITE_FILE),
        "protocol_sha256": sha256(PROTOCOL_FILE),
        "sampling_period": "JJA plus 24-hour buffer",
        "analysis_role": "post-analysis historical extension",
        "source": "Copernicus CDS ERA5-Land hourly time-series point requests",
    })
    encoding = {
        name: {"zlib": True, "complevel": 4, "shuffle": True}
        for name in panel.data_vars
    }
    partial = output.with_suffix(".nc.part")
    partial.unlink(missing_ok=True)
    panel.to_netcdf(partial, engine="h5netcdf", encoding=encoding)
    os.replace(partial, output)
    print(f"[DONE] {output.name}", flush=True)
    return output


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.workers > 12:
        raise SystemExit("--workers must be between 1 and 12")
    if args.max_attempts < 1:
        raise SystemExit("--max-attempts must be positive")
    if not PROTOCOL_FILE.exists():
        raise SystemExit(f"Missing frozen protocol: {PROTOCOL_FILE}")
    sites = pd.read_csv(SITE_FILE).sort_values("site_id").reset_index(drop=True)
    if len(sites) != 121 or sites.site_id.nunique() != 121:
        raise SystemExit("The frozen site manifest is not a 121-site panel")
    selected = sites
    if args.sites:
        selected = sites[sites.site_id.isin(args.sites)].copy()
        missing = sorted(set(args.sites) - set(selected.site_id))
        if missing:
            raise SystemExit(f"Unknown site identifiers: {missing}")

    for directory in (RAW_DIR, SLIM_DIR, OUTPUT_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    print(f"Sites requested: {len(selected)}; workers: {args.workers}", flush=True)
    print(f"Variables: {len(VARIABLES)}; date range: {DATE_RANGE}", flush=True)
    print(f"Protocol SHA-256: {sha256(PROTOCOL_FILE)}", flush=True)

    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_site, row, args.overwrite,
                            args.assemble_only, args.max_attempts): int(row.site_id)
            for _, row in selected.iterrows()
        }
        for future in as_completed(futures):
            site_id = futures[future]
            try:
                future.result()
            except Exception as error:
                failures.append((site_id, repr(error)))
                print(f"[FAIL] site {site_id:03d}: {error}", flush=True)
    if failures:
        raise SystemExit(f"Failed sites: {failures}")
    if len(selected) != 121:
        print("Partial site run complete; yearly panels were not assembled.")
        return

    site_files = [SLIM_DIR / f"site_{site_id:03d}_jja_buffers.nc"
                  for site_id in sites.site_id]
    outputs = [assemble_year(site_files, sites, year, args.overwrite)
               for year in ANALYSIS_YEARS]
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis historical extension",
        "dataset": DATASET,
        "date_range_requested": DATE_RANGE,
        "years_retained": list(ANALYSIS_YEARS),
        "variables": list(VARIABLES),
        "protocol_file": str(PROTOCOL_FILE.relative_to(PROJECT_DIR)),
        "protocol_sha256": sha256(PROTOCOL_FILE),
        "site_manifest": str(SITE_FILE.relative_to(PROJECT_DIR)),
        "site_manifest_sha256": sha256(SITE_FILE),
        "raw_archives_retained": True,
        "archives": [
            {
                "site_id": int(site.site_id),
                "requested_lon": float(site.lon),
                "requested_lat": float(site.lat),
                "file": str((RAW_DIR / f"site_{int(site.site_id):03d}.zip")
                            .relative_to(PROJECT_DIR)),
                "bytes": (RAW_DIR / f"site_{int(site.site_id):03d}.zip").stat().st_size,
                "sha256": sha256(RAW_DIR / f"site_{int(site.site_id):03d}.zip"),
            }
            for site in sites.itertuples(index=False)
        ],
        "outputs": [str(path.relative_to(PROJECT_DIR)) for path in outputs],
    }
    PROVENANCE_FILE.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Provenance: {PROVENANCE_FILE}")


if __name__ == "__main__":
    main()
