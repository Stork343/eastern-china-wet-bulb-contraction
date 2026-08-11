#!/usr/bin/env python3
"""Analyze the frozen 1950--1990 ERA5-Land temporal extension.

The estimator, five physical graph bandwidths, within-record type-7 quartiles,
and product cyclic-shift design follow EXTENSION_ANALYSIS_PROTOCOL.md.  Results
are isolated in output_historical_extension and never alter 1991--2025 output.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_DIR = Path(__file__).resolve().parent.parent
PROTOCOL_FILE = PROJECT_DIR / "EXTENSION_ANALYSIS_PROTOCOL.md"
SITE_FILE = PROJECT_DIR / "data" / "grid" / "eastern_china_121_sites.csv"
DAILY_DIR = PROJECT_DIR / "data" / "era5_historical_extension" / "daily_fields"
OUTPUT_DIR = PROJECT_DIR / "output_historical_extension"
AUDIT_FILE = OUTPUT_DIR / "historical_analysis_manifest.json"
YEARS = tuple(range(1950, 1991))
BANDWIDTHS = np.array([
    125.799765, 251.599530, 503.199060, 1006.398120, 2012.796241
])
RANDOMISATION_DRAWS = 99_999
RANDOMISATION_SEED = 20_260_810


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    lat0 = np.deg2rad(sites.requested_lat.mean())
    coordinates = np.column_stack((
        sites.requested_lon.to_numpy() * 111.32 * np.cos(lat0),
        sites.requested_lat.to_numpy() * 110.57,
    ))
    differences = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(differences).sum(axis=2))


def graph_operators(distance: np.ndarray) -> list[dict]:
    operators = []
    for bandwidth in BANDWIDTHS:
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


def laplacian_metrics(matrix: np.ndarray,
                      operators: list[dict]) -> np.ndarray:
    return np.column_stack([
        np.sum(matrix * (operator["L"] @ matrix), axis=0) /
        (2 * operator["weight_sum"])
        for operator in operators
    ])


def classify_records(day_rows: pd.DataFrame) -> pd.DataFrame:
    out = day_rows.copy()
    # pandas' linear interpolation is Hyndman--Fan type 7, matching R's
    # quantile(..., type=7) default named in the frozen protocol.
    thresholds = out.groupby("record_id").regional_mean_wbt.quantile(
        [0.25, 0.75], interpolation="linear").unstack()
    thresholds.columns = ["q25", "q75"]
    out = out.merge(thresholds, left_on="record_id", right_index=True,
                    validate="many_to_one")
    out["regime"] = np.where(
        out.regional_mean_wbt >= out.q75, "high",
        np.where(out.regional_mean_wbt >= out.q25, "middle", "low"))
    return out


def record_index(day_rows: pd.DataFrame) -> list[dict]:
    records = []
    for record_id, indices in day_rows.groupby(
            "record_id", sort=True).indices.items():
        subset = day_rows.iloc[indices]
        regimes = subset.regime.to_numpy()
        records.append({
            "record_id": int(record_id),
            "year": int(subset.year.iloc[0]),
            "month": int(subset.month.iloc[0]),
            "indices": np.asarray(indices, dtype=int),
            "high": regimes == "high",
            "middle": regimes == "middle",
        })
    return records


def t_summary(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    estimate = float(values.mean())
    standard_error = float(values.std(ddof=1) / math.sqrt(n))
    critical = float(stats.t.ppf(0.975, n - 1)) if n > 1 else math.nan
    statistic = estimate / standard_error if standard_error > 0 else math.nan
    return {
        "years": n,
        "estimate": estimate,
        "standard_error": standard_error,
        "t_statistic": statistic,
        "t_p_value_lower": float(stats.t.cdf(statistic, n - 1)) if n > 1 else math.nan,
        "ci_lower": estimate - critical * standard_error,
        "ci_upper": estimate + critical * standard_error,
        "negative_years": int((values < 0).sum()),
        "negative_fraction": float((values < 0).mean()),
    }


def load_fields() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[Path]]:
    paths = sorted(DAILY_DIR.glob("era5_land_*_jja_daily_fields.csv.gz"))
    if len(paths) != len(YEARS):
        raise RuntimeError(f"Expected 41 historical daily files, found {len(paths)}")
    file_years = tuple(int(path.name.split("_")[2]) for path in paths)
    if file_years != YEARS:
        raise RuntimeError(f"Historical year sequence is not 1950--1990: {file_years}")
    columns = ["year", "month", "record_id", "analysis_date", "day_definition",
               "site_id", "requested_lon", "requested_lat",
               "regional_mean_wbt", "wbt"]
    daily = pd.concat([
        pd.read_csv(path, usecols=columns) for path in paths
    ], ignore_index=True)
    if not daily.day_definition.eq("utc").all():
        raise RuntimeError("Historical extension contains a non-UTC definition")
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = daily[["site_id", "requested_lon", "requested_lat"]].drop_duplicates(
    ).sort_values("site_id").reset_index(drop=True)
    if len(sites) != 121:
        raise RuntimeError("Historical extension does not contain 121 fixed sites")
    manifest = pd.read_csv(SITE_FILE).sort_values("site_id").reset_index(drop=True)
    if not np.array_equal(sites.site_id.to_numpy(), manifest.site_id.to_numpy()):
        raise RuntimeError("Historical site identifiers differ from frozen manifest")
    if not np.allclose(sites[["requested_lon", "requested_lat"]],
                       manifest[["lon", "lat"]], atol=1e-12, rtol=0):
        raise RuntimeError("Historical coordinates differ from frozen manifest")
    counts = daily.groupby(["record_id", "analysis_date"]).site_id.nunique()
    if len(counts) != len(YEARS) * 92 or not counts.eq(121).all():
        raise RuntimeError("Incomplete historical daily spatial fields")
    day_rows = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"))
    day_rows = classify_records(day_rows).sort_values(
        ["record_id", "date"]).reset_index(drop=True)
    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    pivot = pivot.reindex(index=sites.site_id, columns=day_rows.date)
    matrix = pivot.to_numpy()
    if matrix.shape != (121, len(YEARS) * 92) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Unexpected historical matrix shape {matrix.shape}")
    return sites, day_rows, matrix, paths


def record_effects(metrics: np.ndarray, day_rows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for record in record_index(day_rows):
        values = metrics[record["indices"]]
        high = values[record["high"]].mean(axis=0)
        middle = values[record["middle"]].mean(axis=0)
        for j, bandwidth in enumerate(BANDWIDTHS):
            rows.append({
                "record_id": record["record_id"],
                "year": record["year"],
                "month": record["month"],
                "forcing_segment": "1950-1978" if record["year"] <= 1978 else "1979-1990",
                "bandwidth_km": bandwidth,
                "q_middle": middle[j],
                "q_high": high[j],
                "relative_effect": high[j] / middle[j] - 1,
                "n_high": int(record["high"].sum()),
                "n_middle": int(record["middle"].sum()),
            })
    return pd.DataFrame(rows)


def summarize_effects(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    yearly_scale = records.groupby(
        ["year", "forcing_segment", "bandwidth_km"], as_index=False
    ).relative_effect.mean()
    yearly_profile = yearly_scale.groupby(
        ["year", "forcing_segment"], as_index=False
    ).relative_effect.mean().rename(columns={"relative_effect": "profile_effect"})
    if len(yearly_scale) != 41 * 5 or len(yearly_profile) != 41:
        raise RuntimeError("Unequal month/bandwidth aggregation failed")

    overall = t_summary(yearly_profile.profile_effect.to_numpy())
    leave_one_out = np.array([
        yearly_profile.profile_effect.drop(index).mean()
        for index in yearly_profile.index
    ])
    overall.update({
        "analysis_role": "post-analysis historical extension",
        "period": "1950-1990",
        "months_per_summer": 3,
        "bandwidths": 5,
        "loo_min": float(leave_one_out.min()),
        "loo_max": float(leave_one_out.max()),
    })
    pd.DataFrame([overall]).to_csv(
        OUTPUT_DIR / "historical_overall_results.csv", index=False)

    scale_rows = []
    for bandwidth, group in yearly_scale.groupby("bandwidth_km", sort=True):
        scale_rows.append({"bandwidth_km": bandwidth,
                           **t_summary(group.relative_effect.to_numpy())})
    pd.DataFrame(scale_rows).to_csv(
        OUTPUT_DIR / "historical_scale_results.csv", index=False)

    period_rows = []
    for period, group in yearly_profile.groupby("forcing_segment", sort=True):
        period_rows.append({"period": period,
                            **t_summary(group.profile_effect.to_numpy())})
    pd.DataFrame(period_rows).to_csv(
        OUTPUT_DIR / "historical_period_results.csv", index=False)
    records.to_csv(OUTPUT_DIR / "historical_record_scale_effects.csv", index=False)
    yearly_scale.to_csv(
        OUTPUT_DIR / "historical_year_scale_effects.csv", index=False)
    yearly_profile.to_csv(
        OUTPUT_DIR / "historical_year_profile_effects.csv", index=False)
    return yearly_scale, yearly_profile


def product_cyclic_randomisation(metrics: np.ndarray,
                                 day_rows: pd.DataFrame) -> None:
    records = record_index(day_rows)
    if len(records) != 123:
        raise RuntimeError(f"Expected 123 month-year records, found {len(records)}")
    lookups = []
    for record in records:
        values = metrics[record["indices"]]
        n_days = len(values)
        lookup = np.empty((n_days, len(BANDWIDTHS)), dtype=float)
        for offset in range(n_days):
            shifted = values[(np.arange(n_days) + offset) % n_days]
            high = shifted[record["high"]].mean(axis=0)
            middle = shifted[record["middle"]].mean(axis=0)
            lookup[offset] = high / middle - 1
        lookups.append(lookup)
    observed_scale = np.mean([lookup[0] for lookup in lookups], axis=0)
    observed = float(observed_scale.mean())

    rng = np.random.default_rng(RANDOMISATION_SEED)
    null_scale = np.empty((RANDOMISATION_DRAWS, len(BANDWIDTHS)))
    for start in range(0, RANDOMISATION_DRAWS, 5000):
        stop = min(start + 5000, RANDOMISATION_DRAWS)
        total = np.zeros((stop - start, len(BANDWIDTHS)))
        for lookup in lookups:
            offsets = rng.integers(0, len(lookup), size=stop - start)
            total += lookup[offsets]
        null_scale[start:stop] = total / len(lookups)
    null_profile = null_scale.mean(axis=1)
    exceedances = int((null_profile <= observed).sum())
    p_value = (1 + exceedances) / (RANDOMISATION_DRAWS + 1)
    summary = pd.DataFrame([{
        "analysis_role": "post-analysis historical extension",
        "statistic": "equal-record_equal-bandwidth_mean_relative_effect",
        "summers": 41,
        "month_year_records": 123,
        "bandwidths": 5,
        "observed": observed,
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
    null_min = null_scale.min(axis=1)
    for j, bandwidth in enumerate(BANDWIDTHS):
        scale_rows.append({
            "bandwidth_km": bandwidth,
            "observed": observed_scale[j],
            "unadjusted_p": (1 + int((null_scale[:, j] <= observed_scale[j]).sum())) /
            (RANDOMISATION_DRAWS + 1),
            "complete_null_single_step_p": (
                1 + int((null_min <= observed_scale[j]).sum())) /
            (RANDOMISATION_DRAWS + 1),
            "null_mean": null_scale[:, j].mean(),
            "null_sd": null_scale[:, j].std(ddof=1),
        })
    summary.to_csv(OUTPUT_DIR / "historical_global_cyclic_randomisation.csv",
                   index=False)
    pd.DataFrame(scale_rows).to_csv(
        OUTPUT_DIR / "historical_global_cyclic_randomisation_scales.csv", index=False)
    pd.DataFrame({"null_profile_statistic": null_profile}).to_csv(
        OUTPUT_DIR / "historical_global_cyclic_null.csv.gz", index=False,
        compression="gzip")


def energy_decomposition(matrix: np.ndarray, day_rows: pd.DataFrame,
                         operators: list[dict]) -> None:
    months = day_rows.month.to_numpy()
    climatology = np.zeros_like(matrix)
    for month in (6, 7, 8):
        mask = months == month
        climatology[:, mask] = matrix[:, mask].mean(axis=1, keepdims=True)
    anomaly = matrix - climatology
    components = {
        "total": laplacian_metrics(matrix, operators),
        "anomaly_energy": laplacian_metrics(anomaly, operators),
    }
    cross = np.empty_like(components["total"])
    climatology_energy = np.empty_like(components["total"])
    for j, operator in enumerate(operators):
        cross[:, j] = np.sum(
            climatology * (operator["L"] @ anomaly), axis=0
        ) / operator["weight_sum"]
        climatology_energy[:, j] = np.sum(
            climatology * (operator["L"] @ climatology), axis=0
        ) / (2 * operator["weight_sum"])
    components["climatology_anomaly_cross"] = cross
    components["climatology_energy"] = climatology_energy
    identity = np.max(np.abs(
        components["total"] - components["anomaly_energy"] -
        components["climatology_anomaly_cross"] -
        components["climatology_energy"]))
    if identity > 1e-10:
        raise RuntimeError(f"Energy decomposition identity error {identity}")

    rows = []
    for record in record_index(day_rows):
        indices = record["indices"]
        middle_total = components["total"][indices][record["middle"]].mean(axis=0)
        deltas = {
            name: values[indices][record["high"]].mean(axis=0) -
            values[indices][record["middle"]].mean(axis=0)
            for name, values in components.items()
        }
        for j, bandwidth in enumerate(BANDWIDTHS):
            rows.append({
                "year": record["year"], "month": record["month"],
                "record_id": record["record_id"], "bandwidth_km": bandwidth,
                "total_effect": deltas["total"][j] / middle_total[j],
                "anomaly_energy_component": deltas["anomaly_energy"][j] / middle_total[j],
                "climatology_anomaly_cross_component":
                    deltas["climatology_anomaly_cross"][j] / middle_total[j],
                "climatology_energy_component": 0.0,
            })
    records = pd.DataFrame(rows)
    columns = ["total_effect", "anomaly_energy_component",
               "climatology_anomaly_cross_component", "climatology_energy_component"]
    yearly = records.groupby(["year", "bandwidth_km"], as_index=False)[columns].mean()
    summary_rows = []
    for bandwidth, group in yearly.groupby("bandwidth_km", sort=True):
        row = {"bandwidth_km": bandwidth}
        for column in columns:
            for key, value in t_summary(group[column].to_numpy()).items():
                row[f"{column}_{key}"] = value
        row["identity_error"] = abs(
            row["total_effect_estimate"] - row["anomaly_energy_component_estimate"] -
            row["climatology_anomaly_cross_component_estimate"] -
            row["climatology_energy_component_estimate"])
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if summary.identity_error.max() > 1e-10:
        raise RuntimeError("Aggregated historical energy identity failed")
    yearly.to_csv(OUTPUT_DIR / "historical_energy_decomposition_yearly.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "historical_energy_decomposition.csv", index=False)


def basis_decomposition(matrix: np.ndarray, sites: pd.DataFrame,
                        day_rows: pd.DataFrame, operators: list[dict]) -> None:
    latitude = sites.requested_lat.to_numpy()
    longitude = sites.requested_lon.to_numpy()
    bases = {
        "latitude": (latitude - latitude.mean())[:, None],
        "latitude_longitude": np.column_stack((
            latitude - latitude.mean(), longitude - longitude.mean())),
    }
    total = laplacian_metrics(matrix, operators)
    rows = []
    for basis_name, basis in bases.items():
        structured = np.empty_like(total)
        residual = np.empty_like(total)
        for j, operator in enumerate(operators):
            gram = basis.T @ operator["L"] @ basis
            beta = np.linalg.pinv(gram) @ basis.T @ operator["L"] @ matrix
            fitted = basis @ beta
            residual_matrix = matrix - fitted
            structured[:, j] = np.sum(
                fitted * (operator["L"] @ fitted), axis=0
            ) / (2 * operator["weight_sum"])
            residual[:, j] = np.sum(
                residual_matrix * (operator["L"] @ residual_matrix), axis=0
            ) / (2 * operator["weight_sum"])
        identity = np.max(np.abs(total - structured - residual))
        if identity > 1e-10:
            raise RuntimeError(f"{basis_name} historical identity error {identity}")
        for record in record_index(day_rows):
            indices = record["indices"]
            middle_total = total[indices][record["middle"]].mean(axis=0)
            total_delta = total[indices][record["high"]].mean(axis=0) - \
                total[indices][record["middle"]].mean(axis=0)
            structured_delta = structured[indices][record["high"]].mean(axis=0) - \
                structured[indices][record["middle"]].mean(axis=0)
            residual_delta = residual[indices][record["high"]].mean(axis=0) - \
                residual[indices][record["middle"]].mean(axis=0)
            for j, bandwidth in enumerate(BANDWIDTHS):
                rows.append({
                    "basis": basis_name, "year": record["year"],
                    "month": record["month"], "record_id": record["record_id"],
                    "bandwidth_km": bandwidth,
                    "total_effect": total_delta[j] / middle_total[j],
                    "structured_component": structured_delta[j] / middle_total[j],
                    "residual_component": residual_delta[j] / middle_total[j],
                })
    records = pd.DataFrame(rows)
    columns = ["total_effect", "structured_component", "residual_component"]
    yearly = records.groupby(
        ["basis", "year", "bandwidth_km"], as_index=False
    )[columns].mean()
    summary_rows = []
    for (basis, bandwidth), group in yearly.groupby(
            ["basis", "bandwidth_km"], sort=True):
        row = {"basis": basis, "bandwidth_km": bandwidth}
        for column in columns:
            for key, value in t_summary(group[column].to_numpy()).items():
                row[f"{column}_{key}"] = value
        row["identity_error"] = abs(
            row["total_effect_estimate"] - row["structured_component_estimate"] -
            row["residual_component_estimate"])
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if summary.identity_error.max() > 1e-10:
        raise RuntimeError("Aggregated historical basis identity failed")
    yearly.to_csv(OUTPUT_DIR / "historical_basis_decomposition_yearly.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "historical_basis_decomposition.csv", index=False)


def write_manifest(inputs: list[Path], day_rows: pd.DataFrame,
                   matrix: np.ndarray, operators: list[dict]) -> None:
    outputs = sorted(path for path in OUTPUT_DIR.glob("historical_*.csv*")
                     if path != AUDIT_FILE)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis historical extension",
        "protocol_file": str(PROTOCOL_FILE.relative_to(PROJECT_DIR)),
        "protocol_sha256": sha256(PROTOCOL_FILE),
        "site_manifest_sha256": sha256(SITE_FILE),
        "years": [min(YEARS), max(YEARS)],
        "summers": len(YEARS),
        "month_year_records": int(day_rows.record_id.nunique()),
        "daily_fields": int(len(day_rows)),
        "sites": int(matrix.shape[0]),
        "bandwidths_km": BANDWIDTHS.tolist(),
        "quartile_method": "Hyndman-Fan type 7 (linear interpolation)",
        "climatology_reference": {
            "period": "1950-1990",
            "definition": (
                "site-specific mean field by calendar month, estimated only "
                "from the 41 historical summers"
            ),
            "modern_1991_2025_values_used": False,
        },
        "randomisation": {
            "design": "joint product cyclic shift of five dispersion columns within record",
            "draws": RANDOMISATION_DRAWS,
            "seed": RANDOMISATION_SEED,
            "specified_before_historical_value_retrieval": True,
        },
        "forcing_period_check": ["1950-1978", "1979-1990"],
        "forcing_caveat": (
            "ECMWF documents that 1950-1978 ERA5-Land was forced by the "
            "preliminary ERA5 back extension, with sub-optimal representation "
            "of some tropical cyclones; smaller effects are expected over land."
        ),
        "official_documentation": [
            "https://confluence.ecmwf.int/pages/viewpage.action?pageId=501029136",
            "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries",
        ],
        "timeseries_product_note": (
            "The CDS product guide labels this optimized time-series source "
            "experimental and potentially subject to format/source changes."
        ),
        "graph_weight_sums": [operator["weight_sum"] for operator in operators],
        "inputs": [
            {"file": str(path.relative_to(PROJECT_DIR)), "sha256": sha256(path)}
            for path in inputs
        ],
        "outputs": [
            {"file": str(path.relative_to(PROJECT_DIR)), "sha256": sha256(path)}
            for path in outputs
        ],
    }
    AUDIT_FILE.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sites, day_rows, matrix, inputs = load_fields()
    distance = projected_distance(sites)
    operators = graph_operators(distance)
    metrics = laplacian_metrics(matrix, operators)
    records = record_effects(metrics, day_rows)
    summarize_effects(records)
    product_cyclic_randomisation(metrics, day_rows)
    energy_decomposition(matrix, day_rows, operators)
    basis_decomposition(matrix, sites, day_rows, operators)
    write_manifest(inputs, day_rows, matrix, operators)
    print("Historical extension analysis complete.")
    print(pd.read_csv(OUTPUT_DIR / "historical_overall_results.csv").to_string(index=False))
    print(pd.read_csv(OUTPUT_DIR / "historical_period_results.csv").to_string(index=False))
    print(pd.read_csv(
        OUTPUT_DIR / "historical_global_cyclic_randomisation.csv").to_string(index=False))


if __name__ == "__main__":
    main()
