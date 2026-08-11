#!/usr/bin/env python3
"""Verify the nested 465-site spatial-resolution sensitivity bundle."""

from __future__ import annotations

import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import xarray as xr


PROJECT_DIR = Path(__file__).resolve().parent.parent
GRID_DIR = PROJECT_DIR / "data" / "grid"
DENSE_DIR = PROJECT_DIR / "data" / "era5_dense"
PRIMARY_DIR = PROJECT_DIR / "data" / "era5_confirmatory"
OUTPUT_DIR = PROJECT_DIR / "output_dense"
REPORT_FILE = OUTPUT_DIR / "dense_completion_audit.json"
EXPECTED_CORE = {"d2m", "sp", "t2m"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def main() -> None:
    primary = pd.read_csv(GRID_DIR / "eastern_china_121_sites.csv")
    dense = pd.read_csv(GRID_DIR / "eastern_china_dense_sites.csv")
    require(len(primary) == 121, "Primary manifest is not 121 sites")
    require(len(dense) == 465, "Dense manifest is not 465 sites")
    require(int(dense.is_original_site.sum()) == 121,
            "Dense manifest does not mark 121 reused sites")
    nested = dense.loc[dense.is_original_site].merge(
        primary, left_on="original_site_id", right_on="site_id",
        suffixes=("_dense", "_primary"), validate="one_to_one",
    )
    require(
        ((nested.lon_dense - nested.lon_primary).abs() < 1e-12).all() and
        ((nested.lat_dense - nested.lat_primary).abs() < 1e-12).all(),
        "Primary coordinates are not exactly nested",
    )

    new_sites = dense.loc[~dense.is_original_site]
    archives = sorted((DENSE_DIR / "cds_point_archives").glob("site_*.zip"))
    trimmed = sorted((DENSE_DIR / "trimmed_points").glob(
        "site_*_jja_buffers.nc"
    ))
    require(len(archives) == len(new_sites) == 344,
            "Expected 344 new dense archives")
    require(len(trimmed) == 344, "Expected 344 new trimmed point files")
    for path in archives:
        require(zipfile.is_zipfile(path), f"Invalid dense archive: {path.name}")
    dense_by_id = dense.set_index("dense_site_id")
    for path in trimmed:
        dense_id = int(path.name.split("_")[1])
        expected_site = dense_by_id.loc[dense_id]
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 78120,
                    f"Unexpected dense point time size: {path.name}")
            require(set(dataset.data_vars) == EXPECTED_CORE,
                    f"Unexpected dense point variables: {path.name}")
            require(abs(float(dataset.attrs["requested_lon"]) -
                        float(expected_site.lon)) <= 1e-12 and
                    abs(float(dataset.attrs["requested_lat"]) -
                        float(expected_site.lat)) <= 1e-12,
                    f"Dense point coordinate mismatch: {path.name}")

    panels = sorted((DENSE_DIR / "hourly_points").glob(
        "era5_land_*_jja_465sites.nc"
    ))
    require(len(panels) == 35, "Expected 35 dense hourly panels")
    for path in panels:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 2232 and
                    dataset.sizes.get("site") == 465,
                    f"Unexpected dense panel dimensions: {path.name}")
            require(set(dataset.data_vars) == EXPECTED_CORE,
                    f"Unexpected dense panel variables: {path.name}")

    daily = sorted((DENSE_DIR / "daily_fields").glob(
        "era5_land_*_jja_dense_daily_fields.csv.gz"
    ))
    require(len(daily) == 35, "Expected 35 dense daily-field files")
    definitions = {"primary_grid_peak", "dense_grid_peak"}
    for path in daily:
        frame = pd.read_csv(path, usecols=[
            "analysis_definition", "analysis_date", "site_id", "wbt"
        ])
        require(len(frame) == 2 * 92 * 465,
                f"Unexpected dense daily rows: {path.name}")
        require(set(frame.analysis_definition) == definitions,
                f"Unexpected dense definitions: {path.name}")
        counts = frame.groupby("analysis_definition").size()
        require((counts == 92 * 465).all(),
                f"Incomplete dense definition: {path.name}")
        require(frame.wbt.notna().all(), f"Missing dense WBT: {path.name}")

    results = pd.read_csv(OUTPUT_DIR / "dense_primary_results.csv")
    require(set(results.configuration) == {
        "primary_121_reproduction",
        "dense_465_fixed_labels",
        "dense_465_recomputed",
    }, "Dense result configurations are incomplete")
    reproduction = pd.read_csv(
        OUTPUT_DIR / "dense_primary_reproduction_check.csv"
    )
    require(len(reproduction) == 33 and
            reproduction.absolute_error.max() <= 1e-9,
            "Embedded primary-grid reproduction failed")
    identity = pd.read_csv(OUTPUT_DIR / "dense_spatial_identity_check.csv")
    require(len(identity) == 4 and identity.absolute_error.max() <= 1e-10,
            "Dense spatial decomposition identity failed")

    figures = ["fig10_dense_fixed_smooth_surfaces.pdf"]
    for name in figures:
        require((OUTPUT_DIR / name).exists(), f"Missing dense figure: {name}")

    fixed = results.loc[results.configuration == "dense_465_fixed_labels"].iloc[0]
    recomputed = results.loc[
        results.configuration == "dense_465_recomputed"
    ].iloc[0]
    report = {
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_sites": 650,
        "retained_dense_sites": 465,
        "reused_primary_sites": 121,
        "new_archives": len(archives),
        "new_archive_bytes": sum(path.stat().st_size for path in archives),
        "yearly_dense_panels": len(panels),
        "daily_dense_files": len(daily),
        "reproduction_max_error": float(reproduction.absolute_error.max()),
        "spatial_identity_max_error": float(identity.absolute_error.max()),
        "fixed_label_profile_effect": float(fixed.estimate),
        "recomputed_label_profile_effect": float(recomputed.estimate),
        "fixed_label_consistency": bool(fixed.consistency),
        "recomputed_label_consistency": bool(recomputed.consistency),
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
