#!/usr/bin/env python3
"""Post-review spatial, climatology, and station sensitivity analyses.

This script does not replace the protocol-defined primary analysis.  It adds
the targeted checks requested during pre-submission review:

* WGS84 geodesic distances, equal-site versus cosine-latitude weighting,
  inward perturbations of each rectangular boundary, and a non-rectangular
  mainland-China support;
* an exact Eq. 11 climatology--anomaly decomposition using both the original
  inclusive monthly climatology and a leave-one-summer-out climatology; and
* NOAA availability, fixed-support, stricter-day-count, station-defined-event,
  scale-uncertainty, and pressure-conversion checks.

The downloaded dense ERA5-Land domain ends at the original rectangle.  Thus
the boundary analysis can move edges inward but cannot honestly evaluate
outward expansion without new data.  Domain analyses retain the frozen
primary-grid peak hour on each day; both frozen labels and domain-mean labels
are reported so that spatial-support and label changes remain distinguishable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Geod
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
DISCOVERY_YEARS = {2015, 2022}
NOAA_YEARS = (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023, 2025)
H_FACTORS = np.array([0.125, 0.25, 0.5, 1.0, 2.0])
MIN_STATIONS = 10
BROAD_MIN_KM = 500.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-project", type=Path, default=PROJECT,
        help="Project containing the full local data and results directories.",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=PROJECT / "output_revision_sensitivity",
    )
    parser.add_argument(
        "--china-boundary", type=Path,
        help="Optional Natural Earth-compatible country polygon file.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def t_summary(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = len(values)
    if n == 0:
        return {
            "n_years": 0, "estimate": math.nan, "standard_error": math.nan,
            "ci_lower": math.nan, "ci_upper": math.nan,
            "negative_years": 0,
        }
    estimate = float(values.mean())
    if n == 1:
        standard_error = critical = math.nan
    else:
        standard_error = float(values.std(ddof=1) / math.sqrt(n))
        critical = float(stats.t.ppf(0.975, n - 1))
    return {
        "n_years": n,
        "estimate": estimate,
        "standard_error": standard_error,
        "ci_lower": estimate - critical * standard_error if n > 1 else math.nan,
        "ci_upper": estimate + critical * standard_error if n > 1 else math.nan,
        "negative_years": int((values < 0).sum()),
    }


def classify_records(day_rows: pd.DataFrame, mean_column: str) -> pd.DataFrame:
    out = day_rows.copy()
    thresholds = out.groupby("record_id")[mean_column].quantile(
        [0.25, 0.75], interpolation="linear"
    ).unstack()
    thresholds.columns = ["q25", "q75"]
    out = out.drop(columns=["q25", "q75", "regime"], errors="ignore").merge(
        thresholds, left_on="record_id", right_index=True,
        validate="many_to_one",
    )
    out["regime"] = np.where(
        out[mean_column] >= out.q75, "high",
        np.where(out[mean_column] >= out.q25, "middle", "low"),
    )
    return out


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    lat0 = np.deg2rad(sites.requested_lat.mean())
    coordinates = np.column_stack((
        sites.requested_lon.to_numpy() * 111.32 * np.cos(lat0),
        sites.requested_lat.to_numpy() * 110.57,
    ))
    difference = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(difference).sum(axis=2))


def wgs84_distance(sites: pd.DataFrame) -> np.ndarray:
    longitude = sites.requested_lon.to_numpy(dtype=float)
    latitude = sites.requested_lat.to_numpy(dtype=float)
    lon1 = np.broadcast_to(longitude[:, None], (len(sites), len(sites)))
    lat1 = np.broadcast_to(latitude[:, None], (len(sites), len(sites)))
    lon2 = lon1.T
    lat2 = lat1.T
    _, _, metres = Geod(ellps="WGS84").inv(lon1, lat1, lon2, lat2)
    return np.asarray(metres, dtype=float) / 1000.0


def graph_operators(
    distance: np.ndarray,
    bandwidths: np.ndarray,
    node_weights: np.ndarray | None = None,
) -> list[dict]:
    operators = []
    upper = np.triu_indices_from(distance, 1)
    for bandwidth in np.asarray(bandwidths, dtype=float):
        weights = np.exp(-np.square(distance) / (2 * bandwidth**2))
        if node_weights is not None:
            weights *= np.outer(node_weights, node_weights)
        np.fill_diagonal(weights, 0.0)
        edge_weights = weights[upper]
        laplacian = np.diag(weights.sum(axis=1)) - weights
        operators.append({
            "bandwidth_km": float(bandwidth),
            "W": weights,
            "L": laplacian,
            "weight_sum": float(edge_weights.sum()),
            "weighted_mean_pair_distance_km": float(
                np.average(distance[upper], weights=edge_weights)
            ),
            "effective_edges": float(
                edge_weights.sum() ** 2 / np.square(edge_weights).sum()
            ),
        })
    return operators


def field_metrics(matrix: np.ndarray, operators: list[dict]) -> np.ndarray:
    return np.column_stack([
        np.sum(matrix * (operator["L"] @ matrix), axis=0) /
        (2 * operator["weight_sum"])
        for operator in operators
    ])


def record_effects(
    metrics: np.ndarray,
    day_rows: pd.DataFrame,
    bandwidths: np.ndarray,
    analysis: str,
) -> pd.DataFrame:
    rows = []
    for record_id, indices in day_rows.groupby(
        "record_id", sort=True
    ).indices.items():
        subset = day_rows.iloc[indices]
        high = subset.regime.to_numpy() == "high"
        middle = subset.regime.to_numpy() == "middle"
        if high.sum() == 0 or middle.sum() == 0:
            continue
        high_mean = metrics[indices][high].mean(axis=0)
        middle_mean = metrics[indices][middle].mean(axis=0)
        for j, bandwidth in enumerate(bandwidths):
            rows.append({
                "analysis": analysis,
                "record_id": int(record_id),
                "year": int(record_id) // 100,
                "month": int(record_id) % 100,
                "bandwidth_km": float(bandwidth),
                "mean_high": high_mean[j],
                "mean_middle": middle_mean[j],
                "effect": high_mean[j] / middle_mean[j] - 1,
                "n_high": int(high.sum()),
                "n_middle": int(middle.sum()),
            })
    return pd.DataFrame(rows)


def summarize_spatial_curve(
    records: pd.DataFrame,
    operators: list[dict],
    n_sites: int,
    distance_method: str,
    target_weighting: str,
    label_rule: str,
    domain: str,
) -> tuple[pd.DataFrame, dict]:
    yearly = records.groupby(
        ["analysis", "year", "bandwidth_km"], as_index=False
    ).effect.mean()
    held = yearly.loc[~yearly.year.isin(DISCOVERY_YEARS)]
    rows = []
    for j, (bandwidth, group) in enumerate(
        held.groupby("bandwidth_km", sort=True)
    ):
        summary = t_summary(group.effect.to_numpy())
        rows.append({
            "analysis": records.analysis.iloc[0],
            "domain": domain,
            "distance_method": distance_method,
            "target_weighting": target_weighting,
            "label_rule": label_rule,
            "n_sites": n_sites,
            "bandwidth_km": bandwidth,
            "weighted_mean_pair_distance_km":
                operators[j]["weighted_mean_pair_distance_km"],
            "effective_edges": operators[j]["effective_edges"],
            **summary,
        })
    curve = pd.DataFrame(rows)
    profile = held.groupby("year", as_index=False).effect.mean()
    profile_summary = {
        "analysis": records.analysis.iloc[0],
        "domain": domain,
        "distance_method": distance_method,
        "target_weighting": target_weighting,
        "label_rule": label_rule,
        "n_sites": n_sites,
        **t_summary(profile.effect.to_numpy()),
    }
    return curve, profile_summary


def find_results_dir(source: Path) -> Path:
    candidates = [source / "results", source / "output_confirmatory"]
    for candidate in candidates:
        if (candidate / "confirmatory_graph_metadata.csv").exists():
            return candidate
    raise FileNotFoundError("Cannot locate confirmatory result inputs")


def load_primary(
    source: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    directory = source / "data" / "era5_confirmatory" / "daily_fields"
    paths = sorted(directory.glob("era5_land_*_jja_daily_fields.csv.gz"))
    if len(paths) != 35:
        raise RuntimeError(f"Expected 35 primary field files, found {len(paths)}")
    columns = [
        "year", "month", "record_id", "analysis_date", "day_definition",
        "site_id", "requested_lon", "requested_lat", "regional_mean_wbt",
        "wbt",
    ]
    daily = pd.concat([
        pd.read_csv(path, usecols=columns) for path in paths
    ], ignore_index=True)
    daily = daily.loc[daily.day_definition.eq("utc")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = daily[[
        "site_id", "requested_lon", "requested_lat"
    ]].drop_duplicates().sort_values("site_id").reset_index(drop=True)
    day_rows = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"),
    )
    day_rows = classify_records(day_rows, "regional_mean_wbt").sort_values(
        ["record_id", "date"]
    ).reset_index(drop=True)
    matrix = daily.pivot(
        index="site_id", columns="date", values="wbt"
    ).reindex(index=sites.site_id, columns=day_rows.date).to_numpy()
    if matrix.shape != (121, 3220) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid primary matrix {matrix.shape}")
    return sites, day_rows, matrix


def load_dense(
    source: Path, primary_days: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray]:
    directory = source / "data" / "era5_dense" / "daily_fields"
    paths = sorted(directory.glob("era5_land_*_jja_dense_daily_fields.csv.gz"))
    if len(paths) != 35:
        raise RuntimeError(f"Expected 35 dense field files, found {len(paths)}")
    columns = [
        "analysis_date", "analysis_definition", "site_id",
        "requested_lon", "requested_lat", "wbt",
    ]
    daily = pd.concat([
        pd.read_csv(path, usecols=columns) for path in paths
    ], ignore_index=True)
    daily = daily.loc[daily.analysis_definition.eq("primary_grid_peak")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = daily[[
        "site_id", "requested_lon", "requested_lat"
    ]].drop_duplicates().sort_values("site_id").reset_index(drop=True)
    matrix = daily.pivot(
        index="site_id", columns="date", values="wbt"
    ).reindex(index=sites.site_id, columns=primary_days.date).to_numpy()
    if matrix.shape != (465, 3220) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid dense matrix {matrix.shape}")
    return sites, matrix


def locate_china_boundary(explicit: Path | None) -> tuple[Path, str]:
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError(explicit)
        return explicit.resolve(), "user-supplied Natural Earth-compatible file"
    command = [
        "Rscript", "-e",
        "cat(system.file('shapes/world.gpkg', package='spData'))",
    ]
    result = subprocess.run(
        command, check=True, capture_output=True, text=True
    )
    path = Path(result.stdout.strip())
    if not path.exists():
        raise FileNotFoundError(
            "Natural Earth low-resolution world.gpkg from R package spData"
        )
    version = subprocess.run(
        ["Rscript", "-e", "cat(as.character(packageVersion('spData')))"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    return path, f"R spData {version}, Natural Earth low-resolution world.gpkg"


def china_mainland_mask(
    sites: pd.DataFrame, boundary: Path,
) -> np.ndarray:
    import geopandas as gpd

    world = gpd.read_file(boundary)
    if "name_long" in world.columns:
        china = world.loc[world.name_long.eq("China")]
    elif "ADMIN" in world.columns:
        china = world.loc[world.ADMIN.eq("China")]
    else:
        raise RuntimeError("Boundary file has no recognised country-name field")
    if len(china) != 1:
        raise RuntimeError(f"Expected one China feature, found {len(china)}")
    points = gpd.GeoDataFrame(
        sites[["site_id"]].copy(),
        geometry=gpd.points_from_xy(
            sites.requested_lon, sites.requested_lat
        ),
        crs="EPSG:4326",
    )
    china = china.to_crs("EPSG:4326")
    return points.geometry.within(china.geometry.iloc[0]).to_numpy()


def spatial_sensitivity(
    source: Path,
    results_dir: Path,
    boundary_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict], dict]:
    primary_sites, day_rows, primary_matrix = load_primary(source)
    metadata = pd.read_csv(results_dir / "confirmatory_graph_metadata.csv")
    bandwidths = metadata.loc[
        metadata.definition_index.eq(1), ["h_factor", "bandwidth_km"]
    ].drop_duplicates().sort_values("h_factor").bandwidth_km.to_numpy()
    if len(bandwidths) != 5:
        raise RuntimeError("Expected five fixed primary bandwidths")

    curve_parts: list[pd.DataFrame] = []
    profiles: list[dict] = []
    masks: list[dict] = []

    def analyse(
        analysis: str,
        sites: pd.DataFrame,
        matrix: np.ndarray,
        distance: np.ndarray,
        labels: pd.DataFrame,
        domain: str,
        distance_method: str,
        target_weighting: str,
        label_rule: str,
        node_weights: np.ndarray | None = None,
    ) -> None:
        operators = graph_operators(distance, bandwidths, node_weights)
        records = record_effects(
            field_metrics(matrix, operators), labels, bandwidths, analysis
        )
        curve, profile = summarize_spatial_curve(
            records, operators, len(sites), distance_method,
            target_weighting, label_rule, domain,
        )
        curve_parts.append(curve)
        profiles.append(profile)

    projected = projected_distance(primary_sites)
    geodesic = wgs84_distance(primary_sites)
    analyse(
        "primary_equirect_equal_fixed_labels", primary_sites, primary_matrix,
        projected, day_rows, "primary_121", "equirectangular",
        "equal_site", "frozen_primary",
    )
    analyse(
        "primary_wgs84_equal_fixed_labels", primary_sites, primary_matrix,
        geodesic, day_rows, "primary_121", "WGS84_geodesic",
        "equal_site", "frozen_primary",
    )
    area = np.cos(np.deg2rad(primary_sites.requested_lat.to_numpy()))
    analyse(
        "primary_wgs84_area_fixed_labels", primary_sites, primary_matrix,
        geodesic, day_rows, "primary_121", "WGS84_geodesic",
        "cosine_latitude_area", "frozen_primary", area,
    )
    area_days = day_rows.copy()
    area_days["domain_mean_wbt"] = np.average(
        primary_matrix, axis=0, weights=area
    )
    area_days = classify_records(area_days, "domain_mean_wbt")
    analyse(
        "primary_wgs84_area_relabelled", primary_sites, primary_matrix,
        geodesic, area_days, "primary_121", "WGS84_geodesic",
        "cosine_latitude_area", "area_mean_relabelled", area,
    )

    dense_sites, dense_matrix = load_dense(source, day_rows)
    longitude = dense_sites.requested_lon.to_numpy()
    latitude = dense_sites.requested_lat.to_numpy()
    china = china_mainland_mask(dense_sites, boundary_path)
    domain_masks = {
        "full_downloaded_domain": np.ones(len(dense_sites), dtype=bool),
        "south_edge_inward_0.9deg": latitude >= 21.3 - 1e-9,
        "south_edge_inward_1.8deg": latitude >= 22.2 - 1e-9,
        "north_edge_inward_0.9deg": latitude <= 41.1 + 1e-9,
        "north_edge_inward_1.8deg": latitude <= 40.2 + 1e-9,
        "west_edge_inward_0.8deg": longitude >= 105.8 - 1e-9,
        "west_edge_inward_1.6deg": longitude >= 106.6 - 1e-9,
        "east_edge_inward_0.8deg": longitude <= 124.2 + 1e-9,
        "east_edge_inward_1.6deg": longitude <= 123.4 + 1e-9,
        "natural_earth_china_mainland_intersection": china,
    }
    for domain, mask in domain_masks.items():
        for row in dense_sites.assign(in_domain=mask).itertuples(index=False):
            masks.append({
                "domain": domain, "site_id": int(row.site_id),
                "longitude": row.requested_lon,
                "latitude": row.requested_lat,
                "in_domain": bool(row.in_domain),
            })
        sites = dense_sites.loc[mask].reset_index(drop=True)
        matrix = dense_matrix[mask]
        distance = wgs84_distance(sites)
        metrics_operators = graph_operators(distance, bandwidths)
        metrics = field_metrics(matrix, metrics_operators)
        for label_rule, labels in (
            ("frozen_primary", day_rows),
            ("domain_mean_relabelled", None),
        ):
            if labels is None:
                labels = day_rows.copy()
                labels["domain_mean_wbt"] = matrix.mean(axis=0)
                labels = classify_records(labels, "domain_mean_wbt")
            analysis = f"{domain}_wgs84_equal_{label_rule}"
            records = record_effects(
                metrics, labels, bandwidths, analysis
            )
            curve, profile = summarize_spatial_curve(
                records, metrics_operators, len(sites), "WGS84_geodesic",
                "equal_site", label_rule, domain,
            )
            curve_parts.append(curve)
            profiles.append(profile)

        if domain in {
            "full_downloaded_domain",
            "natural_earth_china_mainland_intersection",
        }:
            weights = np.cos(np.deg2rad(sites.requested_lat.to_numpy()))
            labels = day_rows.copy()
            labels["domain_mean_wbt"] = np.average(
                matrix, axis=0, weights=weights
            )
            labels = classify_records(labels, "domain_mean_wbt")
            analyse(
                f"{domain}_wgs84_area_relabelled", sites, matrix, distance,
                labels, domain, "WGS84_geodesic", "cosine_latitude_area",
                "area_mean_relabelled", weights,
            )

    curves = pd.concat(curve_parts, ignore_index=True)
    profile_frame = pd.DataFrame(profiles)
    reference = profile_frame.loc[
        profile_frame.analysis.eq("primary_equirect_equal_fixed_labels"),
        "estimate",
    ].iloc[0]
    expected_primary = pd.read_csv(
        results_dir / "confirmatory_primary_results.csv"
    )
    expected_primary = float(expected_primary.loc[
        expected_primary.day_definition.eq("utc"), "estimate"
    ].iloc[0])
    reproduction_error = abs(reference - expected_primary)
    if reproduction_error > 1e-10:
        raise RuntimeError(
            f"Primary spatial reproduction error {reproduction_error}"
        )
    profile_frame["difference_from_primary_reference"] = (
        profile_frame.estimate - reference
    )
    audit = {
        "fixed_bandwidths_km": bandwidths.tolist(),
        "primary_sites": len(primary_sites),
        "dense_sites": len(dense_sites),
        "china_mainland_sites": int(china.sum()),
        "primary_reproduction_estimate": float(reference),
        "stored_primary_estimate": expected_primary,
        "primary_reproduction_absolute_error": reproduction_error,
        "max_projected_minus_wgs84_distance_km": float(
            np.max(np.abs(projected - geodesic))
        ),
        "boundary_limitation": (
            "Only inward boundary movements are estimable because no ERA5-Land "
            "points were downloaded outside the original rectangle."
        ),
        "time_rule": (
            "All dense-domain checks retain each day's frozen primary-grid "
            "peak hour; fixed and domain-relabelled event definitions are both "
            "reported."
        ),
    }
    return curves, profile_frame, masks, audit


def decomposition_records(
    matrix: np.ndarray,
    day_rows: pd.DataFrame,
    operators: list[dict],
    bandwidths: np.ndarray,
    climatology: np.ndarray,
    analysis: str,
) -> tuple[pd.DataFrame, float, float]:
    anomaly = matrix - climatology
    total = field_metrics(matrix, operators)
    anomaly_energy = field_metrics(anomaly, operators)
    climatology_energy = field_metrics(climatology, operators)
    cross = np.empty_like(total)
    for j, operator in enumerate(operators):
        cross[:, j] = np.sum(
            climatology * (operator["L"] @ anomaly), axis=0
        ) / operator["weight_sum"]
    daily_identity = float(np.max(np.abs(
        total - anomaly_energy - climatology_energy - cross
    )))
    rows = []
    max_climatology_contrast = 0.0
    for record_id, indices in day_rows.groupby(
        "record_id", sort=True
    ).indices.items():
        subset = day_rows.iloc[indices]
        high = subset.regime.to_numpy() == "high"
        middle = subset.regime.to_numpy() == "middle"
        denominator = total[indices][middle].mean(axis=0)
        components = {
            "total_effect": (
                total[indices][high].mean(axis=0) -
                total[indices][middle].mean(axis=0)
            ) / denominator,
            "anomaly_energy_component": (
                anomaly_energy[indices][high].mean(axis=0) -
                anomaly_energy[indices][middle].mean(axis=0)
            ) / denominator,
            "climatology_anomaly_cross_component": (
                cross[indices][high].mean(axis=0) -
                cross[indices][middle].mean(axis=0)
            ) / denominator,
            "climatology_energy_component": (
                climatology_energy[indices][high].mean(axis=0) -
                climatology_energy[indices][middle].mean(axis=0)
            ) / denominator,
        }
        max_climatology_contrast = max(
            max_climatology_contrast,
            float(np.max(np.abs(components["climatology_energy_component"])))
        )
        for j, bandwidth in enumerate(bandwidths):
            rows.append({
                "analysis": analysis,
                "record_id": int(record_id),
                "year": int(record_id) // 100,
                "month": int(record_id) % 100,
                "bandwidth_km": float(bandwidth),
                **{name: value[j] for name, value in components.items()},
            })
    return pd.DataFrame(rows), daily_identity, max_climatology_contrast


def climatology_sensitivity(
    source: Path, results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    sites, day_rows, matrix = load_primary(source)
    metadata = pd.read_csv(results_dir / "confirmatory_graph_metadata.csv")
    bandwidths = metadata.loc[
        metadata.definition_index.eq(1), ["h_factor", "bandwidth_km"]
    ].drop_duplicates().sort_values("h_factor").bandwidth_km.to_numpy()
    operators = graph_operators(projected_distance(sites), bandwidths)
    months = day_rows.month.to_numpy()
    years = day_rows.year.to_numpy()
    inclusive = np.zeros_like(matrix)
    leave_one = np.zeros_like(matrix)
    for month in (6, 7, 8):
        month_mask = months == month
        inclusive[:, month_mask] = matrix[:, month_mask].mean(
            axis=1, keepdims=True
        )
        for year in np.unique(years):
            target = month_mask & (years == year)
            reference = month_mask & (years != year)
            leave_one[:, target] = matrix[:, reference].mean(
                axis=1, keepdims=True
            )
    parts = []
    audit = {}
    for analysis, climatology in (
        ("inclusive_35_summer_monthly_climatology", inclusive),
        ("leave_one_summer_out_monthly_climatology", leave_one),
    ):
        records, daily_error, climatology_contrast = decomposition_records(
            matrix, day_rows, operators, bandwidths, climatology, analysis
        )
        if daily_error > 1e-10 or climatology_contrast > 1e-10:
            raise RuntimeError(
                f"{analysis} decomposition failed its algebraic gate: "
                f"daily={daily_error}, climatology contrast="
                f"{climatology_contrast}"
            )
        parts.append(records)
        audit[analysis] = {
            "maximum_daily_identity_error": daily_error,
            "maximum_record_climatology_energy_contrast": climatology_contrast,
        }
    records = pd.concat(parts, ignore_index=True)
    components = [
        "total_effect", "anomaly_energy_component",
        "climatology_anomaly_cross_component",
        "climatology_energy_component",
    ]
    yearly = records.groupby(
        ["analysis", "year", "bandwidth_km"], as_index=False
    )[components].mean()
    held = yearly.loc[~yearly.year.isin(DISCOVERY_YEARS)]
    curve_rows = []
    for (analysis, bandwidth), group in held.groupby(
        ["analysis", "bandwidth_km"], sort=True
    ):
        for component in components:
            curve_rows.append({
                "analysis": analysis,
                "bandwidth_km": bandwidth,
                "component": component,
                **t_summary(group[component].to_numpy()),
            })
    curve = pd.DataFrame(curve_rows)
    profile_year = held.groupby(
        ["analysis", "year"], as_index=False
    )[components].mean()
    profile_rows = []
    for analysis, group in profile_year.groupby("analysis"):
        for component in components:
            profile_rows.append({
                "analysis": analysis,
                "component": component,
                **t_summary(group[component].to_numpy()),
            })
    return curve, pd.DataFrame(profile_rows), audit


def saturation_vapor_pressure_pa(temperature_k: np.ndarray) -> np.ndarray:
    temperature_c = temperature_k - 273.15
    return 611.2 * np.exp(
        17.67 * temperature_c / (temperature_c + 243.5)
    )


def theta_e(
    temperature_k: np.ndarray,
    dewpoint_k: np.ndarray,
    pressure_pa: np.ndarray,
) -> np.ndarray:
    vapor = saturation_vapor_pressure_pa(dewpoint_k)
    mixing = 0.622 * vapor / (pressure_pa - vapor)
    t_lcl = 1 / (
        1 / (dewpoint_k - 56) + np.log(temperature_k / dewpoint_k) / 800
    ) + 56
    theta_l = (
        temperature_k * (100000 / (pressure_pa - vapor)) ** 0.2854 *
        (temperature_k / t_lcl) ** (0.28 * mixing)
    )
    return theta_l * np.exp(
        (3036 / t_lcl - 1.78) * mixing * (1 + 0.448 * mixing)
    )


def wet_bulb_c(
    temperature_k: np.ndarray,
    dewpoint_k: np.ndarray,
    pressure_pa: np.ndarray,
) -> np.ndarray:
    dewpoint_k = np.minimum(dewpoint_k, temperature_k)
    target = theta_e(temperature_k, dewpoint_k, pressure_pa)
    lower = dewpoint_k.copy()
    upper = temperature_k.copy()
    for _ in range(40):
        midpoint = (lower + upper) / 2
        move = theta_e(midpoint, midpoint, pressure_pa) < target
        lower = np.where(move, midpoint, lower)
        upper = np.where(move, upper, midpoint)
    return (lower + upper) / 2 - 273.15


def load_station_inputs(
    source: Path, results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matched_candidates = [
        results_dir / "noaa_extension_era5_matched.csv.gz",
        source / "output_noaa_extension" / "noaa_extension_era5_matched.csv.gz",
    ]
    matched_path = next((path for path in matched_candidates if path.exists()), None)
    if matched_path is None:
        raise FileNotFoundError("NOAA--ERA matched station panel")
    event_candidates = [
        results_dir / "reviewer_event_manifest.csv",
        results_dir / "sensitivity_event_manifest.csv",
        source / "output_confirmatory" / "sensitivity_event_manifest.csv",
    ]
    event_path = next((path for path in event_candidates if path.exists()), None)
    if event_path is None:
        raise FileNotFoundError("Primary event manifest")
    matched = pd.read_csv(matched_path, dtype={"station": str})
    matched["time_utc"] = pd.to_datetime(matched.time_utc, utc=True).dt.tz_localize(None)
    event = pd.read_csv(event_path)
    event["peak_time"] = pd.to_datetime(event.peak_time, utc=True).dt.tz_localize(None)
    return matched, event


def station_graph_fields(
    matched: pd.DataFrame,
    events: pd.DataFrame,
    bandwidths: np.ndarray,
    value_columns: list[str],
    fixed_support: dict[int, set[str]] | None = None,
) -> pd.DataFrame:
    event_columns = [
        "peak_time", "record_id", "year", "month", "regime"
    ]
    selected = matched.merge(
        events[event_columns], left_on="time_utc", right_on="peak_time",
        how="inner", validate="many_to_one",
    )
    rows = []
    for time, group in selected.groupby("time_utc", sort=True):
        group = group.sort_values("station").drop_duplicates("station")
        year = int(group.year_y.iloc[0])
        if fixed_support is not None:
            support = fixed_support.get(year, set())
            group = group.loc[group.station.isin(support)]
        if len(group) < MIN_STATIONS:
            continue
        sites = group.rename(columns={
            "LONGITUDE": "requested_lon", "LATITUDE": "requested_lat"
        })
        distance = wgs84_distance(sites)
        upper = np.triu_indices(len(group), 1)
        edge_distance = distance[upper]
        for bandwidth in bandwidths:
            weight = np.exp(-np.square(edge_distance) / (2 * bandwidth**2))
            row = {
                "time_utc": time,
                "record_id": int(group.record_id.iloc[0]),
                "year": year,
                "month": int(group.month.iloc[0]),
                "regime": group.regime.iloc[0],
                "bandwidth_km": float(bandwidth),
                "stations": int(len(group)),
            }
            for column in value_columns:
                values = group[column].to_numpy()
                row[f"q_{column}"] = float(
                    np.sum(weight * np.square(
                        values[upper[0]] - values[upper[1]]
                    )) / (2 * weight.sum())
                )
            rows.append(row)
    return pd.DataFrame(rows)


def station_record_effects(
    fields: pd.DataFrame,
    analysis: str,
    value_columns: list[str],
    minimum_days: int,
) -> pd.DataFrame:
    rows = []
    for (record_id, bandwidth), group in fields.groupby(
        ["record_id", "bandwidth_km"], sort=True
    ):
        high = group.regime.eq("high")
        middle = group.regime.eq("middle")
        if high.sum() < minimum_days or middle.sum() < minimum_days:
            continue
        for column in value_columns:
            q_column = f"q_{column}"
            rows.append({
                "analysis": analysis,
                "field": column,
                "minimum_days_per_regime": minimum_days,
                "record_id": int(record_id),
                "year": int(record_id) // 100,
                "month": int(record_id) % 100,
                "bandwidth_km": bandwidth,
                "high_days": int(high.sum()),
                "middle_days": int(middle.sum()),
                "mean_stations": float(group.stations.mean()),
                "effect": group.loc[high, q_column].mean() /
                group.loc[middle, q_column].mean() - 1,
            })
    return pd.DataFrame(rows)


def station_uncertainty(records: pd.DataFrame) -> pd.DataFrame:
    if records.empty:
        return pd.DataFrame()
    yearly = records.groupby(
        ["analysis", "field", "minimum_days_per_regime", "year",
         "bandwidth_km"], as_index=False
    ).agg(
        effect=("effect", "mean"), records=("record_id", "nunique"),
        high_days=("high_days", "sum"), middle_days=("middle_days", "sum"),
        mean_stations=("mean_stations", "mean"),
    )
    rows = []
    for keys, group in yearly.groupby(
        ["analysis", "field", "minimum_days_per_regime", "bandwidth_km"],
        sort=True,
    ):
        analysis, field, minimum_days, bandwidth = keys
        rows.append({
            "analysis": analysis,
            "field": field,
            "minimum_days_per_regime": minimum_days,
            "bandwidth_km": bandwidth,
            "records": int(group.records.sum()),
            "high_days": int(group.high_days.sum()),
            "middle_days": int(group.middle_days.sum()),
            "mean_stations": float(group.mean_stations.mean()),
            **t_summary(group.effect.to_numpy()),
        })
    return pd.DataFrame(rows)


def station_profile_summary(records: pd.DataFrame) -> pd.DataFrame:
    if records.empty:
        return pd.DataFrame()
    broad = records.loc[records.bandwidth_km >= BROAD_MIN_KM]
    yearly = broad.groupby(
        ["analysis", "field", "minimum_days_per_regime", "year"],
        as_index=False,
    ).agg(
        effect=("effect", "mean"), records=("record_id", "nunique"),
        high_days=("high_days", "sum"), middle_days=("middle_days", "sum"),
    )
    rows = []
    for keys, group in yearly.groupby(
        ["analysis", "field", "minimum_days_per_regime"], sort=True
    ):
        analysis, field, minimum_days = keys
        rows.append({
            "analysis": analysis,
            "field": field,
            "minimum_days_per_regime": minimum_days,
            "records": int(group.records.sum()),
            "high_days": int(group.high_days.sum()),
            "middle_days": int(group.middle_days.sum()),
            **t_summary(group.effect.to_numpy()),
        })
    return pd.DataFrame(rows)


def availability_table(
    matched: pd.DataFrame, events: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    scheduled = events.loc[
        events.year.isin(NOAA_YEARS) & events.regime.isin(["high", "middle"]),
        ["peak_time", "record_id", "year", "month", "regime"],
    ].copy()
    counts = matched.groupby("time_utc").station.nunique().rename(
        "station_count"
    )
    scheduled = scheduled.join(counts, on="peak_time")
    scheduled["station_count"] = scheduled.station_count.fillna(0).astype(int)
    scheduled["any_station"] = scheduled.station_count > 0
    scheduled["eligible_ge10"] = scheduled.station_count >= MIN_STATIONS
    rows = []
    groups = [("all_years", 0, scheduled)]
    groups.extend(("year", int(year), group) for year, group in scheduled.groupby("year"))
    for level, year, frame in groups:
        for regime, group in frame.groupby("regime"):
            rows.append({
                "aggregation": level,
                "year": year if year else np.nan,
                "regime": regime,
                "scheduled_days": len(group),
                "days_with_any_station": int(group.any_station.sum()),
                "days_eligible_ge10": int(group.eligible_ge10.sum()),
                "eligibility_fraction_ge10": float(group.eligible_ge10.mean()),
                "station_count_mean": float(group.station_count.mean()),
                "station_count_sd": float(group.station_count.std(ddof=1)),
                "station_count_min": int(group.station_count.min()),
                "station_count_q25": float(group.station_count.quantile(.25)),
                "station_count_median": float(group.station_count.median()),
                "station_count_q75": float(group.station_count.quantile(.75)),
                "station_count_max": int(group.station_count.max()),
            })
    yearly = scheduled.groupby(["year", "regime"]).agg(
        station_count=("station_count", "mean"),
        eligible_fraction=("eligible_ge10", "mean"),
    ).unstack("regime")
    count_difference = (
        yearly[("station_count", "high")] -
        yearly[("station_count", "middle")]
    )
    eligible_difference = (
        yearly[("eligible_fraction", "high")] -
        yearly[("eligible_fraction", "middle")]
    )
    audit = {
        "scheduled_high_middle_fields": len(scheduled),
        "fields_with_any_station": int(scheduled.any_station.sum()),
        "fields_with_at_least_10_stations": int(scheduled.eligible_ge10.sum()),
        "high_minus_middle_station_count": t_summary(count_difference.to_numpy()),
        "high_minus_middle_eligibility_fraction":
            t_summary(eligible_difference.to_numpy()),
    }
    return pd.DataFrame(rows), audit


def fixed_support_by_year(
    matched: pd.DataFrame, events: pd.DataFrame,
) -> tuple[dict[int, set[str]], dict]:
    selected = matched.merge(
        events[["peak_time", "year"]], left_on=["time_utc", "year"],
        right_on=["peak_time", "year"], how="inner",
    )
    field_counts = selected.groupby("time_utc").station.nunique()
    eligible_times = set(field_counts.loc[field_counts >= MIN_STATIONS].index)
    selected = selected.loc[selected.time_utc.isin(eligible_times)]
    supports: dict[int, set[str]] = {}
    details = {}
    global_sets = []
    for year, group in selected.groupby("year"):
        time_sets = [set(x.station) for _, x in group.groupby("time_utc")]
        common = set.intersection(*time_sets) if time_sets else set()
        supports[int(year)] = common
        global_sets.append(common)
        details[str(int(year))] = {
            "eligible_fields": len(time_sets),
            "common_stations": len(common),
            "station_ids": sorted(common),
        }
    global_common = set.intersection(*global_sets) if global_sets else set()
    return supports, {
        "definition": (
            "Within each summer, intersection of station IDs across every "
            "frozen-label event field that originally had at least 10 stations."
        ),
        "year_details": details,
        "global_intersection_of_year_specific_supports": sorted(global_common),
    }


def station_defined_events(matched: pd.DataFrame) -> pd.DataFrame:
    frame = matched.loc[
        matched.year.isin(NOAA_YEARS) & matched.time_utc.dt.month.isin((6, 7, 8))
    ].copy()
    hourly = frame.groupby("time_utc", as_index=False).agg(
        station_count=("station", "nunique"),
        station_mean_wbt=("observed_wbt_c", "mean"),
    )
    hourly = hourly.loc[hourly.station_count >= MIN_STATIONS].copy()
    hourly["date"] = hourly.time_utc.dt.floor("D")
    peaks = hourly.sort_values(
        ["date", "station_mean_wbt", "time_utc"],
        ascending=[True, False, True],
    ).drop_duplicates("date")
    peaks["year"] = peaks.date.dt.year
    peaks["month"] = peaks.date.dt.month
    peaks["record_id"] = peaks.year * 100 + peaks.month
    peaks = peaks.rename(columns={
        "time_utc": "peak_time", "station_mean_wbt": "regional_mean_wbt"
    })
    peaks = classify_records(peaks, "regional_mean_wbt")
    return peaks.loc[peaks.regime.isin(["high", "middle"])]


def station_sensitivity(
    source: Path,
    results_dir: Path,
    bandwidths: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    matched, event = load_station_inputs(source, results_dir)
    matched["observed_wbt_slp_as_surface_c"] = wet_bulb_c(
        matched.temperature_c.to_numpy() + 273.15,
        matched.dewpoint_c.to_numpy() + 273.15,
        matched.slp_hpa.to_numpy() * 100,
    )
    pressure_difference = (
        matched.observed_wbt_slp_as_surface_c - matched.observed_wbt_c
    )
    elevation = matched[["station", "ELEVATION"]].drop_duplicates()
    pressure_audit = {
        "conversion_formula": (
            "p_station_Pa = SLP_hPa * 100 * "
            "max(1 - 2.25577e-5 * elevation_m, 0.1)^5.2559"
        ),
        "station_elevation_min_m": float(elevation.ELEVATION.min()),
        "station_elevation_max_m": float(elevation.ELEVATION.max()),
        "matched_station_hours": len(matched),
        "slp_as_surface_minus_converted_wbt_mean_c":
            float(pressure_difference.mean()),
        "slp_as_surface_minus_converted_wbt_mean_absolute_c":
            float(pressure_difference.abs().mean()),
        "slp_as_surface_minus_converted_wbt_max_absolute_c":
            float(pressure_difference.abs().max()),
    }

    frozen_events = event.loc[
        event.year.isin(NOAA_YEARS) & event.regime.isin(["high", "middle"])
    ].copy()
    availability, availability_audit = availability_table(
        matched, frozen_events
    )
    values = [
        "observed_wbt_c", "observed_wbt_slp_as_surface_c", "era_wbt_c"
    ]
    dynamic_fields = station_graph_fields(
        matched, frozen_events, bandwidths, values
    )
    record_parts = []
    for minimum_days in (1, 2, 3, 5):
        record_parts.append(station_record_effects(
            dynamic_fields,
            f"dynamic_support_frozen_era_labels_min{minimum_days}",
            values, minimum_days,
        ))

    supports, support_audit = fixed_support_by_year(
        matched, frozen_events
    )
    eligible_times = set(dynamic_fields.time_utc.unique())
    common_events = frozen_events.loc[
        frozen_events.peak_time.isin(eligible_times)
    ]
    fixed_fields = station_graph_fields(
        matched, common_events, bandwidths, values, supports
    )
    record_parts.append(station_record_effects(
        fixed_fields, "year_specific_fixed_common_support_frozen_era_labels",
        values, 1,
    ))

    station_events = station_defined_events(matched)
    own_fields = station_graph_fields(
        matched, station_events, bandwidths, values
    )
    record_parts.append(station_record_effects(
        own_fields, "dynamic_support_station_defined_peak_and_labels",
        values, 5,
    ))
    records = pd.concat(record_parts, ignore_index=True)
    uncertainty = station_uncertainty(records)
    profiles = station_profile_summary(records)
    audit = {
        "availability": availability_audit,
        "fixed_support": support_audit,
        "pressure": pressure_audit,
        "station_defined_events": {
            "peak_days": int(station_events.peak_time.nunique()),
            "month_year_records": int(station_events.record_id.nunique()),
            "high_days": int(station_events.regime.eq("high").sum()),
            "middle_days": int(station_events.regime.eq("middle").sum()),
            "daily_peak_rule": (
                "Within each UTC day, maximise the equal-station observed WBT "
                "mean over hours with at least 10 matched stations; apply "
                "record-specific quartile labels to those station means."
            ),
        },
    }
    return uncertainty, profiles, availability, records, audit


def comprehensive_summary(
    spatial_profiles: pd.DataFrame,
    decomposition_profiles: pd.DataFrame,
    station_profiles: pd.DataFrame,
    station_availability: pd.DataFrame,
    station_audit: dict,
) -> pd.DataFrame:
    rows = []
    for row in spatial_profiles.itertuples(index=False):
        rows.append({
            "section": "spatial_support",
            "analysis": row.analysis,
            "field_or_component": "five_scale_raw_ratio",
            "estimate": row.estimate,
            "standard_error": row.standard_error,
            "ci_lower": row.ci_lower,
            "ci_upper": row.ci_upper,
            "n_years": row.n_years,
            "n_sites": row.n_sites,
            "notes": (
                f"{row.domain}; {row.distance_method}; "
                f"{row.target_weighting}; {row.label_rule}"
            ),
        })
    for row in decomposition_profiles.itertuples(index=False):
        rows.append({
            "section": "climatology_decomposition",
            "analysis": row.analysis,
            "field_or_component": row.component,
            "estimate": row.estimate,
            "standard_error": row.standard_error,
            "ci_lower": row.ci_lower,
            "ci_upper": row.ci_upper,
            "n_years": row.n_years,
            "notes": "Five-bandwidth mean; exact Eq. 11 decomposition.",
        })
    for row in station_profiles.itertuples(index=False):
        rows.append({
            "section": "station_comparison",
            "analysis": row.analysis,
            "field_or_component": row.field,
            "estimate": row.estimate,
            "standard_error": row.standard_error,
            "ci_lower": row.ci_lower,
            "ci_upper": row.ci_upper,
            "n_years": row.n_years,
            "n_records": row.records,
            "n_high_days": row.high_days,
            "n_middle_days": row.middle_days,
            "notes": (
                f"Equal mean over bandwidths >=500 km; minimum "
                f"{row.minimum_days_per_regime} retained day(s) per regime."
            ),
        })
    availability_all = station_availability.loc[
        station_availability.aggregation.eq("all_years")
    ]
    for row in availability_all.itertuples(index=False):
        rows.append({
            "section": "station_availability",
            "analysis": "frozen_era_event_schedule",
            "field_or_component": f"{row.regime}_mean_station_count",
            "estimate": row.station_count_mean,
            "n_fields": row.scheduled_days,
            "notes": (
                f"{row.days_eligible_ge10} fields had at least 10 stations; "
                f"eligibility fraction={row.eligibility_fraction_ge10:.6f}."
            ),
        })
    for key, value in station_audit["pressure"].items():
        if isinstance(value, (int, float)):
            rows.append({
                "section": "station_pressure",
                "analysis": "slp_as_surface_vs_standard_atmosphere_conversion",
                "field_or_component": key,
                "estimate": value,
                "notes": station_audit["pressure"]["conversion_formula"],
            })
    return pd.DataFrame(rows)


def methods_note(
    spatial_profiles: pd.DataFrame,
    decomposition_profiles: pd.DataFrame,
    station_profiles: pd.DataFrame,
    station_audit: dict,
    spatial_audit: dict,
    boundary_description: str,
) -> str:
    def spatial_estimate(analysis: str) -> float:
        return 100 * float(spatial_profiles.loc[
            spatial_profiles.analysis.eq(analysis), "estimate"
        ].iloc[0])

    def decomposition_estimate(component: str) -> float:
        return 100 * float(decomposition_profiles.loc[
            decomposition_profiles.analysis.eq(
                "leave_one_summer_out_monthly_climatology"
            ) & decomposition_profiles.component.eq(component),
            "estimate",
        ].iloc[0])

    def station_estimate(analysis: str) -> float:
        return 100 * float(station_profiles.loc[
            station_profiles.analysis.eq(analysis) &
            station_profiles.field.eq("observed_wbt_c"), "estimate"
        ].iloc[0])

    pressure = station_audit["pressure"]
    availability = station_audit["availability"]
    lines = [
        "# Revision sensitivity analysis note",
        "",
        "These are post-review sensitivity analyses; they do not replace the "
        "protocol-defined primary finite-record summary. All reported "
        "intervals are t-scaled between-summer variability intervals, not "
        "confidence intervals under serial dependence.",
        "",
        "## Spatial support and distance",
        "",
        f"The frozen equirectangular/equal-site calculation reproduced the "
        f"stored result to {spatial_audit['primary_reproduction_absolute_error']:.2e}. "
        f"Changing only distance to WGS84 gave {spatial_estimate('primary_wgs84_equal_fixed_labels'):.2f}%, "
        f"versus {spatial_estimate('primary_equirect_equal_fixed_labels'):.2f}% "
        f"for the original calculation. WGS84/cosine-latitude weighting gave "
        f"{spatial_estimate('primary_wgs84_area_fixed_labels'):.2f}% with frozen "
        f"labels and {spatial_estimate('primary_wgs84_area_relabelled'):.2f}% "
        f"after relabelling.",
        "",
        f"The non-rectangular land target uses {boundary_description}: the "
        "country feature named China, excluding the separately stored Taiwan "
        "feature, intersected with the downloaded dense lattice and its frozen "
        f"ERA5-Land-valid mask. It retained {spatial_audit['china_mainland_sites']} "
        f"sites and gave {spatial_estimate('natural_earth_china_mainland_intersection_wgs84_equal_domain_mean_relabelled'):.2f}%. "
        "All one-edge sensitivity analyses move the edge inward; outward "
        "movement is not estimable without new ERA5-Land downloads.",
        "",
        "## Leave-one-summer-out climatology",
        "",
        f"The exact Eq. 11 leave-one-summer-out decomposition gave "
        f"{decomposition_estimate('anomaly_energy_component'):.2f} percentage "
        f"points for anomaly energy and "
        f"{decomposition_estimate('climatology_anomaly_cross_component'):.2f} "
        "percentage points for the climatology--anomaly cross term on the "
        "five-bandwidth average. This is an algebraic structural "
        "interpretation, not a causal mechanism attribution.",
        "",
        "## NOAA station comparison",
        "",
        f"Frozen ERA labels with dynamic support gave "
        f"{station_estimate('dynamic_support_frozen_era_labels_min1'):.2f}% on "
        f"the broad-scale average. A year-specific common-station support gave "
        f"{station_estimate('year_specific_fixed_common_support_frozen_era_labels'):.2f}%, "
        f"and station-defined peak hours and event labels with at least five "
        f"days per regime gave "
        f"{station_estimate('dynamic_support_station_defined_peak_and_labels'):.2f}%.",
        "",
        f"Among {availability['scheduled_high_middle_fields']} frozen high/middle "
        f"event fields, {availability['fields_with_at_least_10_stations']} had at "
        "least 10 stations. The high-minus-middle mean station-count difference "
        f"was {availability['high_minus_middle_station_count']['estimate']:.2f} "
        "stations. The fixed-support analysis uses a separate station "
        "intersection within each summer; the all-summer intersection contains "
        f"only {len(station_audit['fixed_support']['global_intersection_of_year_specific_supports'])} "
        "stations and cannot support a regional graph.",
        "",
        f"Station pressure was computed as `{pressure['conversion_formula']}`. "
        f"Station elevations ranged from {pressure['station_elevation_min_m']:.0f} "
        f"to {pressure['station_elevation_max_m']:.0f} m. Treating sea-level "
        "pressure directly as surface pressure changed WBT by "
        f"{pressure['slp_as_surface_minus_converted_wbt_mean_absolute_c']:.3f} "
        f"degrees C on average in absolute value and at most "
        f"{pressure['slp_as_surface_minus_converted_wbt_max_absolute_c']:.3f} "
        "degrees C.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    source = args.source_project.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results_dir = find_results_dir(source)
    boundary_path, boundary_description = locate_china_boundary(
        args.china_boundary
    )

    spatial_curve, spatial_profiles, domain_masks, spatial_audit = (
        spatial_sensitivity(source, results_dir, boundary_path)
    )
    decomposition_curve, decomposition_profiles, decomposition_audit = (
        climatology_sensitivity(source, results_dir)
    )
    bandwidths = np.sort(spatial_curve.bandwidth_km.unique())
    station_curve, station_profiles, station_availability, station_records, station_audit = (
        station_sensitivity(source, results_dir, bandwidths)
    )
    summary = comprehensive_summary(
        spatial_profiles, decomposition_profiles, station_profiles,
        station_availability, station_audit,
    )

    output_frames = {
        "revision_sensitivity_summary.csv": summary,
        "revision_spatial_scale_curves.csv": spatial_curve,
        "revision_spatial_profile_summary.csv": spatial_profiles,
        "revision_domain_masks.csv": pd.DataFrame(domain_masks),
        "revision_climatology_loso_scale_curve.csv": decomposition_curve,
        "revision_climatology_loso_profile_summary.csv": decomposition_profiles,
        "revision_station_scale_uncertainty.csv": station_curve,
        "revision_station_availability.csv": station_availability,
        "revision_station_record_effects.csv": station_records,
    }
    for name, frame in output_frames.items():
        frame.to_csv(output / name, index=False)
    note_path = output / "REVISION_SENSITIVITY_NOTE.md"
    note_path.write_text(
        methods_note(
            spatial_profiles, decomposition_profiles, station_profiles,
            station_audit, spatial_audit, boundary_description,
        ),
        encoding="utf-8",
    )

    script_path = Path(__file__).resolve()
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": script_path.name,
        "script_sha256": sha256(script_path),
        "analysis_status": "post-review sensitivity; does not redefine primary",
        "source_project": "." if source == PROJECT.resolve() else source.name,
        "output_directory": (
            str(output.relative_to(source)) if output.is_relative_to(source)
            else output.name
        ),
        "spatial": spatial_audit,
        "climatology_decomposition": decomposition_audit,
        "station": station_audit,
        "china_land_operational_definition": {
            "source": boundary_description,
            "file": boundary_path.name,
            "file_sha256": sha256(boundary_path),
            "feature": "country feature named China (ISO CN where available)",
            "taiwan": (
                "Excluded: the source stores Taiwan as a separate feature, and "
                "only the China feature is selected."
            ),
            "intersection": (
                "Point support is the intersection of this polygon, the "
                "downloaded 105--125 E by 20.4--42 N lattice, and the frozen "
                "ERA5-Land-valid mask."
            ),
        },
        "outputs": [],
    }
    for name in output_frames:
        path = output / name
        audit["outputs"].append({
            "file": name, "rows": len(output_frames[name]),
            "sha256": sha256(path),
        })
    audit["outputs"].append({
        "file": note_path.name, "rows": None, "sha256": sha256(note_path),
    })
    audit_path = output / "revision_sensitivity_audit.json"
    audit_path.write_text(
        json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    print("Revision sensitivity analyses complete")
    print(spatial_profiles[[
        "analysis", "n_sites", "estimate", "ci_lower", "ci_upper"
    ]].to_string(index=False))
    print(decomposition_profiles.to_string(index=False))
    print(station_profiles.to_string(index=False))
    print(f"Outputs: {output}")


if __name__ == "__main__":
    main()
