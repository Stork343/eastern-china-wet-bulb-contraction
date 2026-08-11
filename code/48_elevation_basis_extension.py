#!/usr/bin/env python3
"""Acquire invariant ERA5-Land elevation and add it to the spatial basis."""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
PROTOCOL = PROJECT / "EXTENSION_ANALYSIS_PROTOCOL.md"
SITE_FILE = PROJECT / "data" / "grid" / "eastern_china_121_sites.csv"
DAILY_DIR = PROJECT / "data" / "era5_confirmatory" / "daily_fields"
PRIMARY_BASIS = PROJECT / "output_confirmatory" / "extended_basis_decomposition.csv"
BASE = PROJECT / "data" / "era5_invariant"
RAW_FILE = BASE / "era5_land_geopotential_20200101.nc"
ELEVATION_FILE = BASE / "era5_land_121_site_elevation.csv"
OUTPUT = PROJECT / "output_elevation_basis"
RECORD_FILE = OUTPUT / "elevation_basis_record_components.csv"
YEAR_FILE = OUTPUT / "elevation_basis_year_components.csv"
SUMMARY_FILE = OUTPUT / "elevation_basis_summary.csv"
TEX_FILE = OUTPUT / "elevation_basis_table.tex"
AUDIT_FILE = OUTPUT / "elevation_basis_audit.json"
DATASET = "reanalysis-era5-land"
BANDWIDTHS = np.array([125.799765, 251.599530, 503.199060, 1006.398120, 2012.796241])
DISCOVERY = {2015, 2022}
GRAVITY = 9.80665


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire_geopotential() -> None:
    if RAW_FILE.exists() and RAW_FILE.stat().st_size > 1000:
        return
    import cdsapi

    BASE.mkdir(parents=True, exist_ok=True)
    request = {
        "variable": ["geopotential"],
        "year": "2020",
        "month": "01",
        "day": "01",
        "time": "00:00",
        "area": [42, 105, 20, 125],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    partial = RAW_FILE.with_suffix(".download")
    partial.unlink(missing_ok=True)
    cdsapi.Client(quiet=True, progress=False).retrieve(DATASET, request, str(partial))
    if zipfile.is_zipfile(partial):
        with tempfile.TemporaryDirectory(prefix="era5-invariant-") as directory:
            with zipfile.ZipFile(partial) as archive:
                archive.extractall(directory)
            members = sorted(Path(directory).rglob("*.nc"))
            if len(members) != 1:
                raise RuntimeError(f"Expected one NetCDF member, found {len(members)}")
            os.replace(members[0], RAW_FILE)
        partial.unlink(missing_ok=True)
    else:
        os.replace(partial, RAW_FILE)


def extract_elevation() -> pd.DataFrame:
    import xarray as xr

    sites = pd.read_csv(SITE_FILE).sort_values("site_id").reset_index(drop=True)
    with xr.open_dataset(RAW_FILE) as dataset:
        variable = "z" if "z" in dataset.data_vars else "geopotential"
        if variable not in dataset.data_vars:
            raise RuntimeError(f"No geopotential variable in {list(dataset.data_vars)}")
        field = dataset[variable]
        for dimension in list(field.dims):
            if dimension not in {"latitude", "longitude"}:
                field = field.isel({dimension: 0})
        rows = []
        for site in sites.itertuples(index=False):
            value = field.sel(latitude=site.lat, longitude=site.lon, method="nearest")
            matched_lat = float(value.latitude)
            matched_lon = float(value.longitude)
            if abs(matched_lat - site.lat) > 0.051 or abs(matched_lon - site.lon) > 0.051:
                raise RuntimeError(f"Invariant grid mismatch at site {site.site_id}")
            geopotential = float(value)
            rows.append({
                "site_id": int(site.site_id), "lon": float(site.lon), "lat": float(site.lat),
                "matched_lon": matched_lon, "matched_lat": matched_lat,
                "geopotential_m2_s2": geopotential,
                "elevation_m": geopotential / GRAVITY,
            })
    result = pd.DataFrame(rows)
    if len(result) != 121 or not result.elevation_m.between(-100, 6000).all():
        raise RuntimeError("Invalid invariant elevation panel")
    result.to_csv(ELEVATION_FILE, index=False)
    return result


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    lat0 = np.deg2rad(sites.lat.mean())
    coordinates = np.column_stack((sites.lon * 111.32 * np.cos(lat0), sites.lat * 110.57))
    delta = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(delta).sum(axis=2))


def operators(sites: pd.DataFrame) -> list[dict]:
    distance = projected_distance(sites)
    result = []
    for bandwidth in BANDWIDTHS:
        weight = np.exp(-np.square(distance) / (2 * bandwidth**2))
        np.fill_diagonal(weight, 0)
        result.append({
            "bandwidth_km": bandwidth,
            "L": np.diag(weight.sum(axis=1)) - weight,
            "weight_sum": float(np.triu(weight, 1).sum()),
        })
    return result


def metric(matrix: np.ndarray, op: dict) -> np.ndarray:
    return np.sum(matrix * (op["L"] @ matrix), axis=0) / (2 * op["weight_sum"])


def load_daily() -> tuple[pd.DataFrame, np.ndarray]:
    paths = sorted(DAILY_DIR.glob("era5_land_*_jja_daily_fields.csv.gz"))
    columns = ["year", "month", "record_id", "analysis_date", "day_definition",
               "site_id", "regional_mean_wbt", "wbt"]
    daily = pd.concat([pd.read_csv(path, usecols=columns) for path in paths], ignore_index=True)
    daily = daily.loc[daily.day_definition.eq("utc")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    days = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"))
    days = days.sort_values(["record_id", "date"]).reset_index(drop=True)
    quantiles = days.groupby("record_id").regional_mean_wbt.quantile(
        [0.25, 0.75], interpolation="linear").unstack()
    quantiles.columns = ["q25", "q75"]
    days = days.merge(quantiles, left_on="record_id", right_index=True, validate="many_to_one")
    days["regime"] = np.where(
        days.regional_mean_wbt >= days.q75, "high",
        np.where(days.regional_mean_wbt >= days.q25, "middle", "low"))
    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    matrix = pivot.reindex(index=np.arange(1, 122), columns=days.date).to_numpy()
    if matrix.shape != (121, 3220) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid daily matrix {matrix.shape}")
    return days, matrix


def t_summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    n = len(values)
    estimate = float(values.mean())
    se = float(values.std(ddof=1) / math.sqrt(n))
    critical = float(stats.t.ppf(0.975, n - 1))
    return {
        "years": n, "estimate": estimate, "standard_error": se,
        "ci_lower": estimate - critical * se, "ci_upper": estimate + critical * se,
        "negative_years": int((values < 0).sum()),
    }


def analyze(elevation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    days, matrix = load_daily()
    op_list = operators(elevation)
    columns = {
        "latitude": ["lat"],
        "latitude_longitude": ["lat", "lon"],
        "latitude_longitude_elevation": ["lat", "lon", "elevation_m"],
    }
    rows = []
    maximum_identity = 0.0
    for op in op_list:
        total = metric(matrix, op)
        for basis_name, names in columns.items():
            basis = elevation[names].to_numpy(dtype=float)
            basis = (basis - basis.mean(axis=0)) / basis.std(axis=0, ddof=0)
            gram = basis.T @ op["L"] @ basis
            beta = np.linalg.pinv(gram, rcond=1e-12) @ basis.T @ op["L"] @ matrix
            fitted = basis @ beta
            residual_matrix = matrix - fitted
            structured = metric(fitted, op)
            residual = metric(residual_matrix, op)
            maximum_identity = max(maximum_identity, float(np.max(np.abs(total - structured - residual))))
            for record_id, indices in days.groupby("record_id", sort=True).indices.items():
                indices = np.asarray(indices, dtype=int)
                subset = days.iloc[indices]
                high = subset.regime.to_numpy() == "high"
                middle = subset.regime.to_numpy() == "middle"
                denominator = total[indices][middle].mean()
                rows.append({
                    "basis": basis_name, "record_id": int(record_id),
                    "year": int(subset.year.iloc[0]), "month": int(subset.month.iloc[0]),
                    "analysis_role": ("development" if int(subset.year.iloc[0]) in DISCOVERY else "held_out"),
                    "bandwidth_km": op["bandwidth_km"],
                    "total_effect": (total[indices][high].mean() - total[indices][middle].mean()) / denominator,
                    "structured_component": (structured[indices][high].mean() - structured[indices][middle].mean()) / denominator,
                    "residual_component": (residual[indices][high].mean() - residual[indices][middle].mean()) / denominator,
                })
    records = pd.DataFrame(rows)
    years = records.groupby(
        ["basis", "year", "analysis_role", "bandwidth_km"], as_index=False
    )[["total_effect", "structured_component", "residual_component"]].mean()
    held = years.loc[years.analysis_role.eq("held_out")]
    summaries = []
    for (basis, bandwidth), group in held.groupby(["basis", "bandwidth_km"], sort=True):
        row = {"basis": basis, "bandwidth_km": bandwidth}
        for column in ("total_effect", "structured_component", "residual_component"):
            for name, value in t_summary(group[column].to_numpy()).items():
                row[f"{column}_{name}"] = value
        row["identity_error"] = abs(
            row["total_effect_estimate"] - row["structured_component_estimate"] -
            row["residual_component_estimate"])
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    return records, years, summary, maximum_identity


def write_tex(summary: pd.DataFrame) -> None:
    pivot = {}
    for basis, group in summary.groupby("basis"):
        pivot[basis] = group.set_index("bandwidth_km")
    lines = [
        "% Generated by code/48_elevation_basis_extension.py",
        r"\begin{tabular}{rrrrrrrr}", r"\toprule",
        r"Bandwidth & Total & Latitude & Residual & Planar & Residual & Geographic-- & Residual \\",
        r"(km) & & structured & & structured & & topographic & \\", r"\midrule",
    ]
    for bandwidth in BANDWIDTHS:
        total = pivot["latitude"].loc[bandwidth, "total_effect_estimate"]
        values = []
        for basis in ("latitude", "latitude_longitude", "latitude_longitude_elevation"):
            row = pivot[basis].loc[bandwidth]
            values.extend([row.structured_component_estimate, row.residual_component_estimate])
        lines.append(
            f"{bandwidth:,.0f} & {100 * total:.2f} & " +
            " & ".join(f"{100 * value:.2f}" for value in values) + r" \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    TEX_FILE.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    acquire_geopotential()
    elevation = extract_elevation()
    records, years, summary, identity = analyze(elevation)
    records.to_csv(RECORD_FILE, index=False)
    years.to_csv(YEAR_FILE, index=False)
    summary.to_csv(SUMMARY_FILE, index=False)
    write_tex(summary)
    primary = pd.read_csv(PRIMARY_BASIS)
    comparisons = []
    for basis, primary_basis in (("latitude", "latitude"),
                                 ("latitude_longitude", "latitude_longitude")):
        new = summary.loc[summary.basis.eq(basis)].sort_values("bandwidth_km")
        old = primary.loc[primary.basis.eq(primary_basis)].sort_values("bandwidth_km")
        comparisons.append(float(np.max(np.abs(
            new.structured_component_estimate.to_numpy() -
            old.structured_component_estimate.to_numpy()))))
    if max(comparisons) > 1e-10 or identity > 1e-9 or summary.identity_error.max() > 1e-10:
        raise RuntimeError("Elevation-basis reproduction or identity check failed")
    outputs = (ELEVATION_FILE, RECORD_FILE, YEAR_FILE, SUMMARY_FILE, TEX_FILE)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis geographic-topographic basis sensitivity",
        "dataset": DATASET,
        "request": {"variable": "geopotential", "date": "2020-01-01 00:00",
                    "area": [42, 105, 20, 125], "data_format": "netcdf"},
        "raw_file": str(RAW_FILE.relative_to(PROJECT)), "raw_sha256": sha256(RAW_FILE),
        "elevation_conversion": "geopotential / 9.80665 m s^-2",
        "elevation_range_m": [float(elevation.elevation_m.min()), float(elevation.elevation_m.max())],
        "protocol_sha256": sha256(PROTOCOL),
        "basis_spaces": list({"latitude": 1, "latitude_longitude": 2,
                              "latitude_longitude_elevation": 3}),
        "maximum_day_identity_error": identity,
        "maximum_summary_identity_error": float(summary.identity_error.max()),
        "maximum_primary_reproduction_error": max(comparisons),
        "outputs": [{"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
                    for path in outputs],
    }
    AUDIT_FILE.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    topo = summary.loc[summary.basis.eq("latitude_longitude_elevation"), [
        "bandwidth_km", "total_effect_estimate", "structured_component_estimate",
        "residual_component_estimate", "structured_component_negative_years"]]
    print(topo.to_string(index=False))
    print(f"Elevation range: {elevation.elevation_m.min():.1f}--{elevation.elevation_m.max():.1f} m")


if __name__ == "__main__":
    main()
