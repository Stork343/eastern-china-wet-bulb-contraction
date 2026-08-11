#!/usr/bin/env python3
"""Analyses added for the systematic JRSS C manuscript extended analysis.

The script keeps the frozen 121-site, five-bandwidth estimator unchanged and
adds post-analysis diagnostics requested during manuscript review:

1. a product cyclic-shift test over all 99 held-out month-year records;
2. physical-scale and alternative effect-measure summaries;
3. an exact climatology-anomaly graph-energy decomposition;
4. latitude-only and planar Laplacian-basis decompositions;
5. a dense bandwidth profile on the primary grid; and
6. spatial-resolution convergence along two nested refinement paths.

All event labels in these diagnostics are the original UTC labels.  The
spatial-resolution calculation also fixes the primary-grid peak times.  The
new decompositions and dense bandwidth curve are post-analysis summaries and
do not redefine the prespecified estimand.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
DAILY_DIR = PROJECT / "data" / "era5_confirmatory" / "daily_fields"
DENSE_DAILY_DIR = PROJECT / "data" / "era5_dense" / "daily_fields"
GRID_DIR = PROJECT / "data" / "grid"
OUTPUT = PROJECT / "output_confirmatory"
DENSE_OUTPUT = PROJECT / "output_dense"
DISCOVERY_YEARS = {2015, 2022}
H_FACTORS = np.array([0.125, 0.25, 0.5, 1.0, 2.0])
RANDOMISATION_SEED = 20260809
RANDOMISATION_DRAWS = 99_999


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    """Match the equirectangular kilometre projection used in the pipeline."""
    lat0 = np.deg2rad(sites.requested_lat.mean())
    coordinates = np.column_stack((
        sites.requested_lon.to_numpy() * 111.32 * np.cos(lat0),
        sites.requested_lat.to_numpy() * 110.57,
    ))
    differences = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(differences).sum(axis=2))


def graph_operators(distance: np.ndarray,
                    bandwidths: np.ndarray) -> list[dict]:
    operators = []
    for bandwidth in np.asarray(bandwidths, dtype=float):
        weights = np.exp(-np.square(distance) / (2 * bandwidth**2))
        np.fill_diagonal(weights, 0.0)
        laplacian = np.diag(weights.sum(axis=1)) - weights
        operators.append({
            "bandwidth_km": float(bandwidth),
            "W": weights,
            "L": laplacian,
            "weight_sum": float(np.triu(weights, 1).sum()),
        })
    return operators


def laplacian_field_metrics(matrix: np.ndarray,
                            operators: list[dict]) -> np.ndarray:
    """Return day-by-scale graph dispersion for a site-by-day matrix."""
    return np.column_stack([
        np.sum(matrix * (operator["L"] @ matrix), axis=0) /
        (2 * operator["weight_sum"])
        for operator in operators
    ])


def edge_field_metrics(matrix: np.ndarray, distance: np.ndarray,
                       bandwidths: np.ndarray,
                       batch_size: int = 128) -> np.ndarray:
    """Graph dispersion with fixed bandwidths, evaluated in day batches."""
    upper = np.triu_indices_from(distance, 1)
    edge_distance = distance[upper]
    weights = np.exp(
        -np.square(edge_distance[:, None]) /
        (2 * np.square(np.asarray(bandwidths)[None, :]))
    )
    weight_sum = weights.sum(axis=0)
    result = np.empty((matrix.shape[1], len(bandwidths)), dtype=float)
    for start in range(0, matrix.shape[1], batch_size):
        stop = min(start + batch_size, matrix.shape[1])
        difference = (matrix[upper[0], start:stop] -
                      matrix[upper[1], start:stop])
        result[start:stop] = (np.square(difference).T @ weights) / (
            2 * weight_sum)
    return result


def classify_records(day_rows: pd.DataFrame) -> pd.DataFrame:
    out = day_rows.copy()
    thresholds = out.groupby("record_id")["regional_mean_wbt"].quantile(
        [0.25, 0.75], interpolation="linear").unstack()
    thresholds.columns = ["q25", "q75"]
    out = out.merge(thresholds, left_on="record_id", right_index=True,
                    validate="many_to_one")
    out["regime"] = np.where(
        out.regional_mean_wbt >= out.q75, "high",
        np.where(out.regional_mean_wbt >= out.q25, "middle", "low"))
    return out


def t_summary(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    estimate = float(values.mean())
    standard_error = float(values.std(ddof=1) / math.sqrt(n))
    critical = float(stats.t.ppf(0.975, n - 1))
    return {
        "years": n,
        "estimate": estimate,
        "standard_error": standard_error,
        "ci_lower": estimate - critical * standard_error,
        "ci_upper": estimate + critical * standard_error,
        "negative_years": int((values < 0).sum()),
    }


def record_index(day_rows: pd.DataFrame,
                 held_out_only: bool = False) -> list[dict]:
    records = []
    for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
        subset = day_rows.iloc[indices]
        year = int(subset.year.iloc[0])
        if held_out_only and year in DISCOVERY_YEARS:
            continue
        regimes = subset.regime.to_numpy()
        records.append({
            "record_id": int(record_id),
            "year": year,
            "month": int(subset.month.iloc[0]),
            "indices": np.asarray(indices, dtype=int),
            "high": regimes == "high",
            "middle": regimes == "middle",
        })
    return records


def record_metric_table(metrics: np.ndarray, day_rows: pd.DataFrame,
                        bandwidths: np.ndarray,
                        analysis: str) -> pd.DataFrame:
    rows = []
    for record in record_index(day_rows):
        values = metrics[record["indices"]]
        high_mean = values[record["high"]].mean(axis=0)
        middle_mean = values[record["middle"]].mean(axis=0)
        for j, bandwidth in enumerate(bandwidths):
            graph_effect = high_mean[j] / middle_mean[j] - 1
            rows.append({
                "analysis": analysis,
                "record_id": record["record_id"],
                "year": record["year"],
                "month": record["month"],
                "bandwidth_km": float(bandwidth),
                "q_middle": middle_mean[j],
                "q_high": high_mean[j],
                "graph_effect": graph_effect,
                "log_effect": math.log(high_mean[j] / middle_mean[j]),
                "bounded_effect": 2 * (high_mean[j] - middle_mean[j]) /
                (high_mean[j] + middle_mean[j]),
                "rms_middle": math.sqrt(2 * middle_mean[j]),
                "rms_high": math.sqrt(2 * high_mean[j]),
                "rms_effect": math.sqrt(high_mean[j] / middle_mean[j]) - 1,
            })
    return pd.DataFrame(rows)


def year_metric_table(records: pd.DataFrame,
                      columns: list[str]) -> pd.DataFrame:
    return records.groupby(
        ["analysis", "year", "bandwidth_km"], as_index=False
    )[columns].mean()


def load_primary_fields() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    paths = sorted(DAILY_DIR.glob("era5_land_*_jja_daily_fields.csv.gz"))
    if len(paths) != 35:
        raise RuntimeError(f"Expected 35 primary daily files, found {len(paths)}")
    columns = ["year", "month", "record_id", "analysis_date",
               "day_definition", "site_id", "requested_lon",
               "requested_lat", "regional_mean_wbt", "wbt"]
    daily = pd.concat([
        pd.read_csv(path, usecols=columns) for path in paths
    ], ignore_index=True)
    daily = daily.loc[daily.day_definition.eq("utc")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = daily[["site_id", "requested_lon", "requested_lat"]].drop_duplicates(
    ).sort_values("site_id").reset_index(drop=True)
    if len(sites) != 121:
        raise RuntimeError("Primary field does not contain 121 sites")
    day_rows = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"))
    day_rows = classify_records(day_rows).sort_values(
        ["record_id", "date"]).reset_index(drop=True)
    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    pivot = pivot.reindex(index=sites.site_id, columns=day_rows.date)
    matrix = pivot.to_numpy()
    if matrix.shape != (121, 35 * 92) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Unexpected primary matrix shape {matrix.shape}")
    return sites, day_rows, matrix


def product_cyclic_randomisation(metrics: np.ndarray,
                                 day_rows: pd.DataFrame,
                                 bandwidths: np.ndarray) -> None:
    records = record_index(day_rows, held_out_only=True)
    if len(records) != 99:
        raise RuntimeError(f"Expected 99 held-out records, found {len(records)}")

    lookups = []
    for record in records:
        values = metrics[record["indices"]]
        n_days = len(values)
        lookup = np.empty((n_days, values.shape[1]), dtype=float)
        for offset in range(n_days):
            shifted = values[(np.arange(n_days) + offset) % n_days]
            high_mean = shifted[record["high"]].mean(axis=0)
            middle_mean = shifted[record["middle"]].mean(axis=0)
            lookup[offset] = high_mean / middle_mean - 1
        lookups.append(lookup)

    observed_scale = np.mean([lookup[0] for lookup in lookups], axis=0)
    observed_profile = float(observed_scale.mean())
    expected = pd.read_csv(OUTPUT / "confirmatory_primary_results.csv")
    expected = float(expected.loc[expected.day_definition.eq("utc"), "estimate"].iloc[0])
    if abs(observed_profile - expected) > 1e-10:
        raise RuntimeError(
            f"Randomisation statistic {observed_profile} does not reproduce {expected}")

    rng = np.random.default_rng(RANDOMISATION_SEED)
    null_scale = np.empty((RANDOMISATION_DRAWS, len(bandwidths)), dtype=float)
    batch_size = 5_000
    for start in range(0, RANDOMISATION_DRAWS, batch_size):
        stop = min(start + batch_size, RANDOMISATION_DRAWS)
        total = np.zeros((stop - start, len(bandwidths)), dtype=float)
        for lookup in lookups:
            offsets = rng.integers(0, len(lookup), size=stop - start)
            total += lookup[offsets]
        null_scale[start:stop] = total / len(lookups)
    null_profile = null_scale.mean(axis=1)
    exceedances = int((null_profile <= observed_profile).sum())
    p_value = (1 + exceedances) / (RANDOMISATION_DRAWS + 1)
    null_min = null_scale.min(axis=1)

    global_summary = pd.DataFrame([{
        "statistic": "five_scale_finite_record_mean",
        "held_out_summers": 33,
        "month_year_records": len(records),
        "bandwidths": len(bandwidths),
        "observed": observed_profile,
        "alternative": "lower_tail",
        "draws": RANDOMISATION_DRAWS,
        "seed": RANDOMISATION_SEED,
        "random_values_at_or_below_observed": exceedances,
        "p_value": p_value,
        "minimum_attainable_p": 1 / (RANDOMISATION_DRAWS + 1),
        "monte_carlo_se": math.sqrt(p_value * (1 - p_value) /
                                    (RANDOMISATION_DRAWS + 1)),
        "null_mean": null_profile.mean(),
        "null_sd": null_profile.std(ddof=1),
        "null_q001": np.quantile(null_profile, 0.001),
        "null_q01": np.quantile(null_profile, 0.01),
        "null_q05": np.quantile(null_profile, 0.05),
        "null_q50": np.quantile(null_profile, 0.50),
    }])
    scale_rows = []
    for j, bandwidth in enumerate(bandwidths):
        scale_rows.append({
            "bandwidth_km": bandwidth,
            "observed": observed_scale[j],
            "unadjusted_p": (1 + int((null_scale[:, j] <=
                                       observed_scale[j]).sum())) /
            (RANDOMISATION_DRAWS + 1),
            "complete_null_single_step_p": (
                1 + int((null_min <= observed_scale[j]).sum())) /
            (RANDOMISATION_DRAWS + 1),
            "null_mean": null_scale[:, j].mean(),
            "null_sd": null_scale[:, j].std(ddof=1),
        })
    global_summary.to_csv(
        OUTPUT / "extended_global_cyclic_randomisation.csv", index=False)
    pd.DataFrame(scale_rows).to_csv(
        OUTPUT / "extended_global_cyclic_randomisation_scales.csv", index=False)
    pd.DataFrame({"null_profile_statistic": null_profile}).to_csv(
        OUTPUT / "extended_global_cyclic_null.csv.gz", index=False,
        compression="gzip")


def physical_and_effect_measure_summaries(records: pd.DataFrame) -> None:
    columns = ["q_middle", "q_high", "graph_effect", "log_effect",
               "bounded_effect", "rms_middle", "rms_high", "rms_effect"]
    yearly = year_metric_table(records, columns)
    held = yearly.loc[~yearly.year.isin(DISCOVERY_YEARS)].copy()
    scale_rows = []
    for bandwidth, group in held.groupby("bandwidth_km", sort=True):
        graph = t_summary(group.graph_effect.to_numpy())
        q_middle = group.q_middle.mean()
        q_high = group.q_high.mean()
        row = {
            "bandwidth_km": bandwidth,
            "q_middle_c2": q_middle,
            "q_high_c2": q_high,
            "rms_middle_c": math.sqrt(2 * q_middle),
            "rms_high_c": math.sqrt(2 * q_high),
            "mean_record_rms_middle_c": group.rms_middle.mean(),
            "mean_record_rms_high_c": group.rms_high.mean(),
            "graph_effect": graph["estimate"],
            "graph_ci_lower": graph["ci_lower"],
            "graph_ci_upper": graph["ci_upper"],
            "negative_years": graph["negative_years"],
            "rms_effect_mean": group.rms_effect.mean(),
            "rms_effect_from_mean_graph": math.sqrt(1 + graph["estimate"]) - 1,
            "bounded_effect": group.bounded_effect.mean(),
            "log_effect": group.log_effect.mean(),
        }
        scale_rows.append(row)
    pd.DataFrame(scale_rows).to_csv(
        OUTPUT / "extended_scale_physical_summary.csv", index=False)

    profile_year = held.groupby("year", as_index=False)[
        ["graph_effect", "log_effect", "bounded_effect", "rms_effect"]
    ].mean()
    rows = []
    for measure in ["graph_effect", "log_effect", "bounded_effect", "rms_effect"]:
        summary = t_summary(profile_year[measure].to_numpy())
        if measure == "log_effect":
            summary.update(
                transformed_estimate=math.exp(summary["estimate"]) - 1,
                transformed_ci_lower=math.exp(summary["ci_lower"]) - 1,
                transformed_ci_upper=math.exp(summary["ci_upper"]) - 1,
            )
        rows.append({"measure": measure, **summary})
    pd.DataFrame(rows).to_csv(
        OUTPUT / "extended_effect_measure_summary.csv", index=False)


def energy_decomposition(matrix: np.ndarray, day_rows: pd.DataFrame,
                         operators: list[dict],
                         bandwidths: np.ndarray) -> None:
    months = day_rows.month.to_numpy()
    climatology = np.zeros_like(matrix)
    for month in (6, 7, 8):
        mask = months == month
        climatology[:, mask] = matrix[:, mask].mean(axis=1, keepdims=True)
    anomaly = matrix - climatology

    component_metrics: dict[str, np.ndarray] = {
        "total": laplacian_field_metrics(matrix, operators),
        "anomaly_energy": laplacian_field_metrics(anomaly, operators),
    }
    cross = np.empty_like(component_metrics["total"])
    climatology_energy = np.empty_like(component_metrics["total"])
    for j, operator in enumerate(operators):
        cross[:, j] = np.sum(
            climatology * (operator["L"] @ anomaly), axis=0
        ) / operator["weight_sum"]
        climatology_energy[:, j] = np.sum(
            climatology * (operator["L"] @ climatology), axis=0
        ) / (2 * operator["weight_sum"])
    component_metrics["climatology_anomaly_cross"] = cross
    component_metrics["climatology_energy"] = climatology_energy
    identity = np.max(np.abs(
        component_metrics["total"] - component_metrics["anomaly_energy"] -
        component_metrics["climatology_anomaly_cross"] -
        component_metrics["climatology_energy"]
    ))
    if identity > 1e-10:
        raise RuntimeError(f"Energy decomposition identity error {identity}")

    rows = []
    for record in record_index(day_rows):
        indices = record["indices"]
        total = component_metrics["total"][indices]
        middle_total = total[record["middle"]].mean(axis=0)
        deltas = {
            component: values[indices][record["high"]].mean(axis=0) -
            values[indices][record["middle"]].mean(axis=0)
            for component, values in component_metrics.items()
        }
        for j, bandwidth in enumerate(bandwidths):
            rows.append({
                "record_id": record["record_id"],
                "year": record["year"],
                "month": record["month"],
                "bandwidth_km": bandwidth,
                "total_effect": deltas["total"][j] / middle_total[j],
                "anomaly_energy_component":
                    deltas["anomaly_energy"][j] / middle_total[j],
                "climatology_anomaly_cross_component":
                    deltas["climatology_anomaly_cross"][j] / middle_total[j],
                # The monthly climatology is identical on high and middle
                # days within a record, so this contrast is zero by
                # construction.  Write the mathematical zero rather than
                # retaining subtraction noise at roughly 1e-16.
                "climatology_energy_component": 0.0,
            })
    records = pd.DataFrame(rows)
    components = ["total_effect", "anomaly_energy_component",
                  "climatology_anomaly_cross_component",
                  "climatology_energy_component"]
    yearly = records.groupby(
        ["year", "bandwidth_km"], as_index=False
    )[components].mean()
    yearly["analysis_role"] = np.where(
        yearly.year.isin(DISCOVERY_YEARS), "development", "held_out")
    held = yearly.loc[yearly.analysis_role.eq("held_out")]
    summary_rows = []
    for bandwidth, group in held.groupby("bandwidth_km", sort=True):
        row = {"bandwidth_km": bandwidth}
        for component in components:
            summary = t_summary(group[component].to_numpy())
            for key, value in summary.items():
                row[f"{component}_{key}"] = value
        row["identity_error"] = abs(
            row["total_effect_estimate"] -
            row["anomaly_energy_component_estimate"] -
            row["climatology_anomaly_cross_component_estimate"] -
            row["climatology_energy_component_estimate"])
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if summary.identity_error.max() > 1e-10:
        raise RuntimeError("Aggregated energy decomposition does not close")
    yearly.to_csv(OUTPUT / "extended_energy_decomposition_yearly.csv", index=False)
    summary.to_csv(OUTPUT / "extended_energy_decomposition.csv", index=False)


def basis_decomposition(matrix: np.ndarray, sites: pd.DataFrame,
                        day_rows: pd.DataFrame, operators: list[dict],
                        bandwidths: np.ndarray) -> None:
    latitude = sites.requested_lat.to_numpy()
    longitude = sites.requested_lon.to_numpy()
    bases = {
        "latitude": (latitude - latitude.mean())[:, None],
        "latitude_longitude": np.column_stack((
            latitude - latitude.mean(), longitude - longitude.mean())),
    }
    total_metrics = laplacian_field_metrics(matrix, operators)
    rows = []
    for basis_name, basis in bases.items():
        structured = np.empty_like(total_metrics)
        residual = np.empty_like(total_metrics)
        max_identity_error = 0.0
        for j, operator in enumerate(operators):
            gram = basis.T @ operator["L"] @ basis
            beta = np.linalg.pinv(gram) @ basis.T @ operator["L"] @ matrix
            fitted = basis @ beta
            structured[:, j] = np.sum(
                fitted * (operator["L"] @ fitted), axis=0
            ) / (2 * operator["weight_sum"])
            residual_matrix = matrix - fitted
            residual[:, j] = np.sum(
                residual_matrix * (operator["L"] @ residual_matrix), axis=0
            ) / (2 * operator["weight_sum"])
            max_identity_error = max(
                max_identity_error,
                float(np.max(np.abs(total_metrics[:, j] -
                                    structured[:, j] - residual[:, j]))))
        if max_identity_error > 1e-10:
            raise RuntimeError(
                f"{basis_name} projection identity error {max_identity_error}")

        for record in record_index(day_rows):
            indices = record["indices"]
            middle_total = total_metrics[indices][record["middle"]].mean(axis=0)
            total_delta = (total_metrics[indices][record["high"]].mean(axis=0) -
                           total_metrics[indices][record["middle"]].mean(axis=0))
            structured_delta = (
                structured[indices][record["high"]].mean(axis=0) -
                structured[indices][record["middle"]].mean(axis=0))
            residual_delta = (
                residual[indices][record["high"]].mean(axis=0) -
                residual[indices][record["middle"]].mean(axis=0))
            for j, bandwidth in enumerate(bandwidths):
                rows.append({
                    "basis": basis_name,
                    "record_id": record["record_id"],
                    "year": record["year"],
                    "month": record["month"],
                    "bandwidth_km": bandwidth,
                    "total_effect": total_delta[j] / middle_total[j],
                    "structured_component": structured_delta[j] / middle_total[j],
                    "residual_component": residual_delta[j] / middle_total[j],
                })
    records = pd.DataFrame(rows)
    basis_components = ["total_effect", "structured_component",
                        "residual_component"]
    yearly = records.groupby(
        ["basis", "year", "bandwidth_km"], as_index=False
    )[basis_components].mean()
    yearly["analysis_role"] = np.where(
        yearly.year.isin(DISCOVERY_YEARS), "development", "held_out")
    held = yearly.loc[yearly.analysis_role.eq("held_out")]
    summary_rows = []
    for (basis, bandwidth), group in held.groupby(
            ["basis", "bandwidth_km"], sort=True):
        row = {"basis": basis, "bandwidth_km": bandwidth}
        for component in ["total_effect", "structured_component",
                          "residual_component"]:
            summary = t_summary(group[component].to_numpy())
            for key, value in summary.items():
                row[f"{component}_{key}"] = value
        row["identity_error"] = abs(
            row["total_effect_estimate"] -
            row["structured_component_estimate"] -
            row["residual_component_estimate"])
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if summary.identity_error.max() > 1e-10:
        raise RuntimeError("Aggregated basis decomposition does not close")
    yearly.to_csv(OUTPUT / "extended_basis_decomposition_yearly.csv", index=False)
    summary.to_csv(OUTPUT / "extended_basis_decomposition.csv", index=False)


def dense_bandwidth_profile(matrix: np.ndarray, distance: np.ndarray,
                            day_rows: pd.DataFrame,
                            fixed_bandwidths: np.ndarray) -> None:
    dense_bandwidths = np.exp(np.linspace(
        np.log(fixed_bandwidths.min()), np.log(fixed_bandwidths.max()), 31))
    dense_metrics = edge_field_metrics(matrix, distance, dense_bandwidths)
    dense_records = record_metric_table(
        dense_metrics, day_rows, dense_bandwidths, "dense_bandwidth_curve")
    dense_yearly = year_metric_table(dense_records, ["graph_effect"])
    dense_held = dense_yearly.loc[~dense_yearly.year.isin(DISCOVERY_YEARS)]
    rows = []
    for bandwidth, group in dense_held.groupby("bandwidth_km", sort=True):
        rows.append({"source": "31_point_curve", "bandwidth_km": bandwidth,
                     **t_summary(group.graph_effect.to_numpy())})

    fixed = pd.read_csv(OUTPUT / "confirmatory_scale_results.csv")
    metadata = pd.read_csv(OUTPUT / "confirmatory_graph_metadata.csv")
    metadata = metadata.loc[metadata.definition_index.eq(1),
                            ["metric", "bandwidth_km"]].drop_duplicates()
    fixed = fixed.loc[
        fixed.day_definition.eq("utc") & fixed.metric.str.startswith("graph_h_")
    ].merge(metadata, on="metric", validate="one_to_one")
    for row in fixed.itertuples(index=False):
        rows.append({
            "source": "prespecified_scale",
            "bandwidth_km": row.bandwidth_km,
            "years": row.years,
            "estimate": row.estimate,
            "standard_error": row.standard_error,
            "ci_lower": row.ci_lower,
            "ci_upper": row.ci_upper,
            "negative_years": row.negative_years,
        })
    pd.DataFrame(rows).to_csv(
        OUTPUT / "extended_dense_bandwidth_profile.csv", index=False)


def load_dense_primary_peak(day_rows: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    paths = sorted(DENSE_DAILY_DIR.glob(
        "era5_land_*_jja_dense_daily_fields.csv.gz"))
    if len(paths) != 35:
        raise RuntimeError(f"Expected 35 dense daily files, found {len(paths)}")
    columns = ["analysis_date", "analysis_definition", "site_id",
               "requested_lon", "requested_lat", "wbt"]
    daily = pd.concat([
        pd.read_csv(path, usecols=columns) for path in paths
    ], ignore_index=True)
    daily = daily.loc[daily.analysis_definition.eq("primary_grid_peak")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = daily[["site_id", "requested_lon", "requested_lat"]].drop_duplicates(
    ).sort_values("site_id").reset_index(drop=True)
    if len(sites) != 465:
        raise RuntimeError(f"Expected 465 dense sites, found {len(sites)}")
    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    pivot = pivot.reindex(index=sites.site_id, columns=day_rows.date)
    matrix = pivot.to_numpy()
    if matrix.shape != (465, 35 * 92) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Unexpected dense matrix shape {matrix.shape}")
    return sites, matrix


def spatial_convergence(day_rows: pd.DataFrame,
                        fixed_bandwidths: np.ndarray) -> None:
    dense_sites, dense_matrix = load_dense_primary_peak(day_rows)
    primary_manifest = pd.read_csv(GRID_DIR / "eastern_china_121_sites.csv")
    dense_manifest = pd.read_csv(GRID_DIR / "eastern_china_dense_sites.csv")
    dense_manifest = dense_manifest.sort_values("dense_site_id").reset_index(drop=True)
    if not np.array_equal(dense_manifest.dense_site_id.to_numpy(),
                          dense_sites.site_id.to_numpy()):
        raise RuntimeError("Dense site manifest and daily fields disagree")

    primary_lons = primary_manifest.lon.unique()
    primary_lats = primary_manifest.lat.unique()
    lon_mask = np.any(np.isclose(
        dense_manifest.lon.to_numpy()[:, None], primary_lons[None, :]), axis=1)
    lat_mask = np.any(np.isclose(
        dense_manifest.lat.to_numpy()[:, None], primary_lats[None, :]), axis=1)
    primary_mask = dense_manifest.is_original_site.astype(bool).to_numpy()
    configurations = {
        "primary_121": primary_mask,
        "longitude_refined": lat_mask,
        "latitude_refined": lon_mask,
        "dense_465": np.ones(len(dense_manifest), dtype=bool),
    }
    if not all(np.all(mask[primary_mask]) for mask in configurations.values()):
        raise RuntimeError("A convergence grid does not contain every primary site")

    primary_scale = pd.read_csv(OUTPUT / "confirmatory_scale_results.csv")
    primary_meta = pd.read_csv(OUTPUT / "confirmatory_graph_metadata.csv")
    primary_meta = primary_meta.loc[primary_meta.definition_index.eq(1),
                                    ["metric", "bandwidth_km"]].drop_duplicates()
    primary_scale = primary_scale.loc[
        primary_scale.day_definition.eq("utc") &
        primary_scale.metric.str.startswith("graph_h_")
    ].merge(primary_meta, on="metric", validate="one_to_one")
    dense_scale = pd.read_csv(DENSE_OUTPUT / "dense_scale_results.csv")
    dense_scale = dense_scale.loc[
        dense_scale.configuration.eq("dense_465_fixed_labels"),
        ["bandwidth_km", "years", "estimate", "standard_error", "ci_lower",
         "ci_upper", "negative_years"]
    ]

    rows = []
    for name, mask in configurations.items():
        n_sites = int(mask.sum())
        if name == "primary_121":
            endpoint = primary_scale
            for row in endpoint.itertuples(index=False):
                rows.append({
                    "configuration": name, "sites": n_sites,
                    "bandwidth_km": row.bandwidth_km,
                    "years": row.years, "estimate": row.estimate,
                    "standard_error": row.standard_error,
                    "ci_lower": row.ci_lower,
                    "ci_upper": row.ci_upper,
                    "negative_years": row.negative_years,
                })
            continue
        if name == "dense_465":
            for row in dense_scale.itertuples(index=False):
                rows.append({
                    "configuration": name, "sites": n_sites,
                    "bandwidth_km": row.bandwidth_km,
                    "years": row.years, "estimate": row.estimate,
                    "standard_error": row.standard_error,
                    "ci_lower": row.ci_lower,
                    "ci_upper": row.ci_upper,
                    "negative_years": row.negative_years,
                })
            continue

        sites = dense_sites.loc[mask].reset_index(drop=True)
        matrix = dense_matrix[mask]
        distance = projected_distance(sites)
        metrics = edge_field_metrics(matrix, distance, fixed_bandwidths)
        records = record_metric_table(
            metrics, day_rows, fixed_bandwidths, name)
        yearly = year_metric_table(records, ["graph_effect"])
        held = yearly.loc[~yearly.year.isin(DISCOVERY_YEARS)]
        for bandwidth, group in held.groupby("bandwidth_km", sort=True):
            rows.append({"configuration": name, "sites": n_sites,
                         "bandwidth_km": bandwidth,
                         **t_summary(group.graph_effect.to_numpy())})

    convergence = pd.DataFrame(rows)
    # CSV round-tripping changes the last few binary digits of the endpoint
    # bandwidths, whereas the intermediate grids use the in-memory values.
    # A six-decimal physical-distance key joins the same fixed bandwidths
    # without treating numerically identical scales as different.
    convergence["bandwidth_key"] = convergence.bandwidth_km.round(6)
    reference = convergence.loc[convergence.configuration.eq("dense_465"),
                                ["bandwidth_key", "estimate"]].rename(
                                    columns={"estimate": "dense_reference"})
    convergence = convergence.merge(reference, on="bandwidth_key",
                                    validate="many_to_one")
    convergence["difference_from_dense"] = (
        convergence.estimate - convergence.dense_reference)
    summary = convergence.groupby(
        ["configuration", "sites"], as_index=False
    ).difference_from_dense.agg(
        max_absolute_difference=lambda x: np.abs(x).max(),
        rmse=lambda x: math.sqrt(np.square(x).mean()))
    convergence.to_csv(
        DENSE_OUTPUT / "extended_spatial_convergence.csv", index=False)
    summary.to_csv(
        DENSE_OUTPUT / "extended_spatial_convergence_summary.csv", index=False)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    DENSE_OUTPUT.mkdir(parents=True, exist_ok=True)
    sites, day_rows, matrix = load_primary_fields()
    distance = projected_distance(sites)
    median_distance = np.median(distance[np.tril_indices_from(distance, -1)])
    fixed_bandwidths = H_FACTORS * median_distance
    operators = graph_operators(distance, fixed_bandwidths)
    metrics = laplacian_field_metrics(matrix, operators)
    records = record_metric_table(
        metrics, day_rows, fixed_bandwidths, "primary_raw_wbt")

    product_cyclic_randomisation(metrics, day_rows, fixed_bandwidths)
    physical_and_effect_measure_summaries(records)
    energy_decomposition(matrix, day_rows, operators, fixed_bandwidths)
    basis_decomposition(matrix, sites, day_rows, operators, fixed_bandwidths)
    dense_bandwidth_profile(matrix, distance, day_rows, fixed_bandwidths)
    spatial_convergence(day_rows, fixed_bandwidths)

    audit = {
        "script": Path(__file__).name,
        "primary_sites": int(len(sites)),
        "days": int(matrix.shape[1]),
        "held_out_summers": 33,
        "held_out_records": 99,
        "fixed_bandwidths_km": fixed_bandwidths.tolist(),
        "randomisation_draws": RANDOMISATION_DRAWS,
        "randomisation_seed": RANDOMISATION_SEED,
        "outputs": [
            "extended_global_cyclic_randomisation.csv",
            "extended_global_cyclic_randomisation_scales.csv",
            "extended_global_cyclic_null.csv.gz",
            "extended_scale_physical_summary.csv",
            "extended_effect_measure_summary.csv",
            "extended_energy_decomposition.csv",
            "extended_energy_decomposition_yearly.csv",
            "extended_basis_decomposition.csv",
            "extended_basis_decomposition_yearly.csv",
            "extended_dense_bandwidth_profile.csv",
            "../output_dense/extended_spatial_convergence.csv",
            "../output_dense/extended_spatial_convergence_summary.csv",
        ],
    }
    (OUTPUT / "extended_analysis_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print("extended analysis analyses complete")
    print(pd.read_csv(
        OUTPUT / "extended_global_cyclic_randomisation.csv").to_string(index=False))
    print(pd.read_csv(
        OUTPUT / "extended_energy_decomposition.csv").to_string(index=False))
    print(pd.read_csv(
        DENSE_OUTPUT / "extended_spatial_convergence_summary.csv").to_string(
            index=False))


if __name__ == "__main__":
    main()
