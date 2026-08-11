#!/usr/bin/env python3
"""Acquire the new sites for the 465-site nested ERA5-Land grid.

The 121 primary locations reuse their immutable confirmatory archives and
trimmed point files. CDS requests are submitted only for the 344 new sites.
Dense outputs live under data/era5_dense and cannot overwrite the primary
confirmatory bundle.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
CODE_DIR = Path(__file__).resolve().parent
SITE_FILE = PROJECT_DIR / "data" / "grid" / "eastern_china_dense_sites.csv"
BASE_DIR = PROJECT_DIR / "data" / "era5_dense"
RAW_DIR = BASE_DIR / "cds_point_archives"
SLIM_DIR = BASE_DIR / "trimmed_points"
OUTPUT_DIR = BASE_DIR / "hourly_points"
PROVENANCE_FILE = BASE_DIR / "cds_dense_provenance.json"

PRIMARY_BASE = PROJECT_DIR / "data" / "era5_confirmatory"
PRIMARY_RAW_DIR = PRIMARY_BASE / "cds_point_archives"
PRIMARY_SLIM_DIR = PRIMARY_BASE / "trimmed_points"
EXPECTED_SITES = 465
CORE_NAMES = ("d2m", "sp", "t2m")


def load_primary_downloader():
    path = CODE_DIR / "21_download_confirmatory_cds_points.py"
    spec = importlib.util.spec_from_file_location("primary_point_downloader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.RAW_DIR = RAW_DIR
    module.SLIM_DIR = SLIM_DIR
    module.OUTPUT_DIR = OUTPUT_DIR
    return module


PRIMARY = load_primary_downloader()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Concurrent CDS requests; keep at 1 unless the account limit changes.",
    )
    parser.add_argument("--max-attempts", type=int, default=12)
    parser.add_argument(
        "--sites", type=int, nargs="+",
        help="Dense site IDs to process. Primary sites in the list are reused.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--assemble-only", action="store_true")
    parser.add_argument("--no-assemble", action="store_true")
    return parser.parse_args()


def dense_row(row: pd.Series) -> pd.Series:
    out = row.copy()
    out["site_id"] = int(row.dense_site_id)
    return out


def trimmed_path(row) -> Path:
    if bool(row.is_original_site):
        return PRIMARY_SLIM_DIR / (
            f"site_{int(row.original_site_id):03d}_jja_buffers.nc"
        )
    return SLIM_DIR / f"site_{int(row.dense_site_id):03d}_jja_buffers.nc"


def archive_path(row) -> Path:
    if bool(row.is_original_site):
        return PRIMARY_RAW_DIR / f"site_{int(row.original_site_id):03d}.zip"
    return RAW_DIR / f"site_{int(row.dense_site_id):03d}.zip"


def validate_reused_site(row) -> Path:
    path = trimmed_path(row)
    archive = archive_path(row)
    if not path.exists() or not PRIMARY.valid_zip(archive):
        raise FileNotFoundError(
            f"Primary cache missing for dense site {int(row.dense_site_id)}"
        )
    return path


def process_new_site(row: pd.Series, args: argparse.Namespace) -> Path:
    return PRIMARY.process_site(
        dense_row(row), PRIMARY.CORE_VARIABLES, args.overwrite,
        args.assemble_only, args.max_attempts,
    )


def assemble_year(site_files: list[Path], sites: pd.DataFrame, year: int,
                  overwrite: bool) -> Path:
    import xarray as xr

    output = OUTPUT_DIR / f"era5_land_{year}_jja_{EXPECTED_SITES}sites.nc"
    if output.exists() and not overwrite:
        print(f"[SKIP] {output.name}", flush=True)
        return output
    start = f"{year}-05-31 16:00:00"
    end = f"{year}-09-01 15:00:00"
    pieces = []
    for path, site in zip(site_files, sites.itertuples(index=False)):
        with xr.open_dataset(path, engine="h5netcdf") as source:
            piece = source[list(CORE_NAMES)].sel(time=slice(start, end)).load()
        piece = piece.expand_dims(site=[int(site.dense_site_id)])
        pieces.append(piece)
    panel = xr.concat(pieces, dim="site", coords="minimal", compat="override")
    panel = panel.assign_coords(
        site_id=("site", sites.dense_site_id.to_numpy(dtype=np.int32)),
        requested_lon=("site", sites.lon.to_numpy()),
        requested_lat=("site", sites.lat.to_numpy()),
        original_site_id=(
            "site", sites.original_site_id.fillna(-1).to_numpy(dtype=np.int32)
        ),
    )
    if panel.sizes.get("time") != 2232 or \
            panel.sizes.get("site") != EXPECTED_SITES:
        raise ValueError(f"Unexpected dimensions for {year}: {dict(panel.sizes)}")
    panel.attrs.update({
        "spatial_manifest": SITE_FILE.name,
        "spatial_manifest_sha256": PRIMARY.sha256(SITE_FILE),
        "sampling_period": "JJA plus 24-hour buffer for complete UTC+8 days",
        "analysis_role": (
            "discovery" if year in PRIMARY.DISCOVERY_YEARS else "confirmatory"
        ),
        "analysis_status": "secondary nested-grid resolution sensitivity",
        "source": "Copernicus CDS ERA5-Land time-series point requests",
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


def write_provenance(sites: pd.DataFrame, outputs: list[Path]) -> None:
    entries = []
    for site in sites.itertuples(index=False):
        archive = archive_path(site)
        entries.append({
            "dense_site_id": int(site.dense_site_id),
            "original_site_id": (
                int(site.original_site_id)
                if not pd.isna(site.original_site_id) else None
            ),
            "reused_primary_archive": bool(site.is_original_site),
            "requested_lon": float(site.lon),
            "requested_lat": float(site.lat),
            "file": str(archive.relative_to(PROJECT_DIR)),
            "bytes": archive.stat().st_size,
            "sha256": PRIMARY.sha256(archive),
        })
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_status": "secondary nested-grid resolution sensitivity",
        "dataset": PRIMARY.DATASET,
        "date_range_requested": PRIMARY.DATE_RANGE,
        "years_retained": list(PRIMARY.ANALYSIS_YEARS),
        "variables": list(PRIMARY.CORE_VARIABLES),
        "site_manifest": str(SITE_FILE.relative_to(PROJECT_DIR)),
        "site_manifest_sha256": PRIMARY.sha256(SITE_FILE),
        "dense_sites": EXPECTED_SITES,
        "reused_primary_sites": int(sites.is_original_site.sum()),
        "new_sites": int((~sites.is_original_site).sum()),
        "archives": entries,
        "outputs": [str(path.relative_to(PROJECT_DIR)) for path in outputs],
    }
    PROVENANCE_FILE.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Provenance: {PROVENANCE_FILE}", flush=True)


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.workers > 12:
        raise SystemExit("--workers must be between 1 and 12")
    if args.max_attempts < 1:
        raise SystemExit("--max-attempts must be positive")

    sites = pd.read_csv(SITE_FILE).sort_values("dense_site_id").reset_index(drop=True)
    sites["is_original_site"] = sites.is_original_site.astype(bool)
    if len(sites) != EXPECTED_SITES or sites.dense_site_id.nunique() != EXPECTED_SITES:
        raise SystemExit(f"Dense manifest is not {EXPECTED_SITES} sites")
    if int(sites.is_original_site.sum()) != 121:
        raise SystemExit("Dense manifest does not contain 121 primary sites")

    selected = sites
    if args.sites:
        selected = sites[sites.dense_site_id.isin(args.sites)].copy()
        missing = sorted(set(args.sites) - set(selected.dense_site_id))
        if missing:
            raise SystemExit(f"Unknown dense site identifiers: {missing}")

    for directory in (RAW_DIR, SLIM_DIR, OUTPUT_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    reused = selected[selected.is_original_site]
    for row in reused.itertuples(index=False):
        validate_reused_site(row)
    new_sites = selected[~selected.is_original_site]
    print(
        f"Selected dense sites: {len(selected)}; reused: {len(reused)}; "
        f"new requests: {len(new_sites)}; workers: {args.workers}",
        flush=True,
    )
    print(
        f"Variables: {len(PRIMARY.CORE_VARIABLES)}; date range: {PRIMARY.DATE_RANGE}",
        flush=True,
    )

    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_new_site, row, args): int(row.dense_site_id)
            for _, row in new_sites.iterrows()
        }
        for future in as_completed(futures):
            site_id = futures[future]
            try:
                future.result()
            except Exception as error:
                failures.append((site_id, repr(error)))
                print(f"[FAIL] dense site {site_id:03d}: {error}", flush=True)
    if failures:
        raise SystemExit(f"Failed dense sites: {failures}")

    if args.no_assemble or len(selected) != EXPECTED_SITES:
        print("Partial dense-site run complete; yearly panels not assembled.")
        return
    site_files = [trimmed_path(row) for row in sites.itertuples(index=False)]
    missing_files = [path for path in site_files if not path.exists()]
    if missing_files:
        raise SystemExit(f"Missing {len(missing_files)} trimmed dense-site files")
    outputs = [
        assemble_year(site_files, sites, year, args.overwrite)
        for year in PRIMARY.ANALYSIS_YEARS
    ]
    write_provenance(sites, outputs)


if __name__ == "__main__":
    main()
