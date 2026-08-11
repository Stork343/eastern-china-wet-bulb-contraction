#!/usr/bin/env python3
"""sensitivity diagnostics for the held-out humid-heat analysis.

This script does not redefine the prespecified estimator.  It calculates the
post-analysis checks needed to delimit its interpretation: log-ratio and
denominator behaviour, calendar-time dependence, spatial weighting and kernel
sensitivity, fixed-hour fields, anomaly fields, a Laplacian latitude-gradient
decomposition, event reproducibility, and station-level agreement.
"""

from __future__ import annotations

import importlib.util
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
DAILY_DIR = PROJECT / "data" / "era5_confirmatory" / "daily_fields"
HOURLY_DIR = PROJECT / "data" / "era5_confirmatory" / "hourly_points"
OUTPUT = PROJECT / "output_confirmatory"
DISCOVERY = {2015, 2022}
H_FACTORS = np.array([0.125, 0.25, 0.5, 1.0, 2.0])


def load_wbt_function():
    path = PROJECT / "code" / "22_build_confirmatory_fields.py"
    spec = importlib.util.spec_from_file_location("field_builder", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.bolton_wbt_c


BOLTON_WBT = load_wbt_function()


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    lat0 = np.deg2rad(sites.requested_lat.mean())
    coordinates = np.column_stack((
        sites.requested_lon.to_numpy() * 111.32 * np.cos(lat0),
        sites.requested_lat.to_numpy() * 110.57,
    ))
    differences = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(differences).sum(axis=2))


def graph_operator(weights: np.ndarray) -> dict[str, np.ndarray | float]:
    weights = np.asarray(weights, dtype=float).copy()
    np.fill_diagonal(weights, 0.0)
    laplacian = np.diag(weights.sum(axis=1)) - weights
    return {
        "W": weights,
        "L": laplacian,
        "weight_sum": float(np.triu(weights, 1).sum()),
    }


def gaussian_operators(distance: np.ndarray) -> list[dict]:
    median_distance = np.median(distance[np.tril_indices_from(distance, -1)])
    out = []
    for factor in H_FACTORS:
        bandwidth = factor * median_distance
        operator = graph_operator(np.exp(-np.square(distance) /
                                         (2 * bandwidth**2)))
        operator.update(h_factor=factor, bandwidth_km=bandwidth,
                        kernel="Gaussian")
        out.append(operator)
    return out


def weighted_quantile(values: np.ndarray, weights: np.ndarray,
                      probability: float) -> float:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights) / weights.sum()
    return float(np.interp(probability, cumulative, values))


def operator_metadata(operators: list[dict], distance: np.ndarray,
                      n_sites: int) -> pd.DataFrame:
    upper = np.triu_indices_from(distance, 1)
    d = distance[upper]
    rows = []
    nearest = np.min(np.where(distance > 0, distance, np.inf), axis=1)
    for operator in operators:
        w = operator["W"][upper]
        degree = operator["W"].sum(axis=1)
        node_effective = np.square(degree) / np.square(operator["W"]).sum(axis=1)
        rows.append({
            "kernel": operator.get("kernel", "unknown"),
            "h_factor": operator["h_factor"],
            "bandwidth_km": operator["bandwidth_km"],
            "nearest_neighbour_km_median": np.median(nearest),
            "weighted_distance_mean_km": np.average(d, weights=w),
            "weighted_distance_median_km": weighted_quantile(d, w, 0.5),
            "weighted_distance_q10_km": weighted_quantile(d, w, 0.1),
            "weighted_distance_q90_km": weighted_quantile(d, w, 0.9),
            "effective_edges": np.square(w.sum()) / np.square(w).sum(),
            "mean_effective_neighbours": np.mean(node_effective),
            "total_edges": n_sites * (n_sites - 1) / 2,
        })
    return pd.DataFrame(rows)


def field_metrics(matrix: np.ndarray, operators: list[dict]) -> np.ndarray:
    """Return day-by-scale graph dispersion for site-by-day fields."""
    return np.column_stack([
        np.sum(matrix * (operator["L"] @ matrix), axis=0) /
        (2 * operator["weight_sum"])
        for operator in operators
    ])


def classify_records(day_table: pd.DataFrame) -> pd.DataFrame:
    out = day_table.copy()
    thresholds = out.groupby("record_id")["regional_mean_wbt"].quantile(
        [0.25, 0.75], interpolation="linear").unstack()
    thresholds.columns = ["q25", "q75"]
    out = out.drop(columns=["q25", "q75"], errors="ignore").merge(
        thresholds, left_on="record_id", right_index=True, validate="many_to_one")
    out["regime"] = np.where(
        out.regional_mean_wbt >= out.q75, "high",
        np.where(out.regional_mean_wbt >= out.q25, "middle", "low"))
    return out


def record_contrasts(metrics: np.ndarray, day_rows: pd.DataFrame,
                     labels: list[str], analysis: str) -> pd.DataFrame:
    rows = []
    for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
        subset = day_rows.iloc[indices]
        high = subset.regime.to_numpy() == "high"
        middle = subset.regime.to_numpy() == "middle"
        for j, label in enumerate(labels):
            high_mean = metrics[indices, j][high].mean()
            middle_mean = metrics[indices, j][middle].mean()
            rows.append({
                "analysis": analysis,
                "record_id": int(record_id),
                "year": int(record_id) // 100,
                "month": int(record_id) % 100,
                "scale": label,
                "mean_high": high_mean,
                "mean_middle": middle_mean,
                "ratio_effect": high_mean / middle_mean - 1,
                "log_effect": np.log(high_mean) - np.log(middle_mean),
                "n_high": int(high.sum()),
                "n_middle": int(middle.sum()),
            })
    return pd.DataFrame(rows)


def year_effects(records: pd.DataFrame, column: str) -> pd.DataFrame:
    scale_years = records.groupby(
        ["analysis", "year", "scale"], as_index=False)[column].mean()
    return scale_years.groupby(["analysis", "year"], as_index=False)[column].mean()


def t_summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    estimate = values.mean()
    se = values.std(ddof=1) / math.sqrt(n)
    critical = stats.t.ppf(0.975, n - 1)
    return {
        "n_years": n,
        "estimate": estimate,
        "se": se,
        "ci_lower": estimate - critical * se,
        "ci_upper": estimate + critical * se,
        "p_one_sided": stats.t.cdf(estimate / se, n - 1),
    }


def summarize_analysis(records: pd.DataFrame, column: str) -> pd.DataFrame:
    yearly = year_effects(records, column)
    rows = []
    for analysis, group in yearly.groupby("analysis"):
        held = group.loc[~group.year.isin(DISCOVERY), column].to_numpy()
        row = {"analysis": analysis, "estimand": column, **t_summary(held)}
        if column == "log_effect":
            row.update(
                transformed_effect=np.exp(row["estimate"]) - 1,
                transformed_ci_lower=np.exp(row["ci_lower"]) - 1,
                transformed_ci_upper=np.exp(row["ci_upper"]) - 1,
            )
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_scales(records: pd.DataFrame) -> pd.DataFrame:
    yearly = records.groupby(
        ["analysis", "year", "scale"], as_index=False).ratio_effect.mean()
    rows = []
    for (analysis, scale), group in yearly.loc[
            ~yearly.year.isin(DISCOVERY)].groupby(["analysis", "scale"]):
        rows.append({"analysis": analysis, "scale": scale,
                     **t_summary(group.ratio_effect.to_numpy())})
    return pd.DataFrame(rows)


def matched_kernel(distance: np.ndarray, target_mean: float,
                   kernel: str) -> tuple[np.ndarray, float]:
    upper = np.triu_indices_from(distance, 1)
    d = distance[upper]

    def weights(parameter: float) -> np.ndarray:
        if kernel == "Exponential":
            return np.exp(-distance / parameter)
        return np.square(np.maximum(1 - distance / parameter, 0))

    lower = 1e-3 if kernel == "Exponential" else np.min(d) * 1.000001
    upper_parameter = np.max(d) * 100
    for _ in range(100):
        midpoint = math.sqrt(lower * upper_parameter)
        w = weights(midpoint)[upper]
        mean_distance = np.average(d, weights=w)
        if mean_distance < target_mean:
            lower = midpoint
        else:
            upper_parameter = midpoint
    parameter = math.sqrt(lower * upper_parameter)
    return weights(parameter), parameter


def calendar_hac(years: np.ndarray, values: np.ndarray, lag: int = 2) -> dict:
    order = np.argsort(years)
    years = years[order].astype(int)
    values = values[order]
    centred = values - values.mean()
    gamma0 = np.mean(np.square(centred))
    long_run = gamma0
    covariances = []
    pairs = []
    value_by_year = dict(zip(years, centred))
    for k in range(1, lag + 1):
        products = [value_by_year[y] * value_by_year[y + k]
                    for y in years if y + k in value_by_year]
        gamma = float(np.mean(products))
        long_run += 2 * (1 - k / (lag + 1)) * gamma
        covariances.append(gamma)
        pairs.append(len(products))
    se = math.sqrt(max(long_run, np.finfo(float).eps) / len(values))
    critical = stats.t.ppf(0.975, len(values) - 1)
    estimate = values.mean()
    return {
        "estimate": estimate,
        "se": se,
        "ci_lower": estimate - critical * se,
        "ci_upper": estimate + critical * se,
        "p_one_sided": stats.t.cdf(estimate / se, len(values) - 1),
        "gamma0": gamma0,
        "gamma1": covariances[0],
        "gamma2": covariances[1],
        "pairs_lag1": pairs[0],
        "pairs_lag2": pairs[1],
    }


def compressed_hac(values: np.ndarray, lag: int = 2) -> dict:
    years = np.arange(len(values))
    return calendar_hac(years, values, lag)


def temporal_diagnostics(yearly: pd.DataFrame) -> pd.DataFrame:
    held = yearly.loc[~yearly.year.isin(DISCOVERY)].sort_values("year")
    years = held.year.to_numpy()
    values = held.ratio_effect.to_numpy()
    rows = []
    for name, summary in (("compressed_lag", compressed_hac(values)),
                          ("calendar_lag", calendar_hac(years, values))):
        rows.append({"diagnostic": name, **summary})
    centred_year = years - years.mean()
    trend = stats.linregress(centred_year, values)
    rows.append({
        "diagnostic": "linear_trend",
        "estimate": trend.slope,
        "se": trend.stderr,
        "ci_lower": trend.slope - stats.t.ppf(0.975, len(values) - 2) * trend.stderr,
        "ci_upper": trend.slope + stats.t.ppf(0.975, len(values) - 2) * trend.stderr,
        "p_two_sided": trend.pvalue,
    })
    early = values[years <= 2007]
    late = values[years >= 2008]
    difference = late.mean() - early.mean()
    pooled_se = math.sqrt(early.var(ddof=1) / len(early) +
                          late.var(ddof=1) / len(late))
    welch_df = (early.var(ddof=1) / len(early) + late.var(ddof=1) / len(late))**2 / (
        (early.var(ddof=1) / len(early))**2 / (len(early) - 1) +
        (late.var(ddof=1) / len(late))**2 / (len(late) - 1))
    critical = stats.t.ppf(0.975, welch_df)
    rows.append({
        "diagnostic": "late_minus_early",
        "estimate": difference,
        "se": pooled_se,
        "ci_lower": difference - critical * pooled_se,
        "ci_upper": difference + critical * pooled_se,
        "p_two_sided": 2 * stats.t.sf(abs(difference / pooled_se), welch_df),
        "early_mean": early.mean(),
        "late_mean": late.mean(),
    })
    negative = int(np.sum(values < 0))
    proportion = negative / len(values)
    z = stats.norm.ppf(0.975)
    denominator = 1 + z**2 / len(values)
    centre = (proportion + z**2 / (2 * len(values))) / denominator
    half = z * math.sqrt(proportion * (1 - proportion) / len(values) +
                         z**2 / (4 * len(values)**2)) / denominator
    rows.append({
        "diagnostic": "negative_year_fraction",
        "estimate": proportion,
        "ci_lower": centre - half,
        "ci_upper": centre + half,
        "negative_years": negative,
        "p_one_sided_independent_sign_reference": stats.binom.sf(negative - 1,
                                                                  len(values), 0.5),
    })
    return pd.DataFrame(rows)


def gradient_decomposition(
        matrix: np.ndarray, day_rows: pd.DataFrame,
        operators: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    latitude = day_rows.attrs["sites"].requested_lat.to_numpy()
    ell = latitude - latitude.mean()
    components = []
    for operator in operators:
        denominator = float(ell @ operator["L"] @ ell)
        beta = (ell @ operator["L"] @ matrix) / denominator
        gradient = np.square(beta) * denominator / (2 * operator["weight_sum"])
        total = np.sum(matrix * (operator["L"] @ matrix), axis=0) / (
            2 * operator["weight_sum"])
        components.append((total, gradient, total - gradient))
    rows = []
    for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
        subset = day_rows.iloc[indices]
        high = subset.regime.to_numpy() == "high"
        middle = subset.regime.to_numpy() == "middle"
        for j, (total, gradient, residual) in enumerate(components):
            baseline = total[indices][middle].mean()
            total_effect = (total[indices][high].mean() -
                            total[indices][middle].mean()) / baseline
            gradient_effect = (gradient[indices][high].mean() -
                               gradient[indices][middle].mean()) / baseline
            residual_effect = (residual[indices][high].mean() -
                               residual[indices][middle].mean()) / baseline
            rows.append({
                "record_id": int(record_id),
                "year": int(record_id) // 100,
                "month": int(record_id) % 100,
                "bandwidth_km": operators[j]["bandwidth_km"],
                "total_effect": total_effect,
                "gradient_component": gradient_effect,
                "residual_component": residual_effect,
            })
    record = pd.DataFrame(rows)
    year = record.groupby(["year", "bandwidth_km"], as_index=False).mean(numeric_only=True)
    held = year.loc[~year.year.isin(DISCOVERY)]
    summaries = []
    for bandwidth, group in held.groupby("bandwidth_km", sort=True):
        total = t_summary(group.total_effect.to_numpy())
        gradient = t_summary(group.gradient_component.to_numpy())
        residual = t_summary(group.residual_component.to_numpy())
        summaries.append({
            "bandwidth_km": bandwidth,
            "total_effect": total["estimate"],
            "total_ci_lower": total["ci_lower"],
            "total_ci_upper": total["ci_upper"],
            "total_negative_years": int((group.total_effect < 0).sum()),
            "gradient_component": gradient["estimate"],
            "gradient_ci_lower": gradient["ci_lower"],
            "gradient_ci_upper": gradient["ci_upper"],
            "gradient_negative_years": int(
                (group.gradient_component < 0).sum()),
            "residual_component": residual["estimate"],
            "residual_ci_lower": residual["ci_lower"],
            "residual_ci_upper": residual["ci_upper"],
            "residual_positive_years": int(
                (group.residual_component > 0).sum()),
        })
    summary = pd.DataFrame(summaries)
    summary["gradient_share"] = (summary.gradient_component /
                                 summary.total_effect)
    summary["identity_error"] = np.abs(
        summary.total_effect - summary.gradient_component -
        summary.residual_component)
    year["analysis_role"] = np.where(
        year.year.isin(DISCOVERY), "development", "held_out")
    return summary, year


def anomaly_matrices(matrix: np.ndarray, day_rows: pd.DataFrame,
                     sites: pd.DataFrame) -> dict[str, np.ndarray]:
    months = day_rows.month.to_numpy()
    years = day_rows.year.to_numpy()
    climatology = np.zeros_like(matrix)
    standardised = np.zeros_like(matrix)
    within_record = np.zeros_like(matrix)
    for month in (6, 7, 8):
        mask = months == month
        means = matrix[:, mask].mean(axis=1, keepdims=True)
        standard_deviation = matrix[:, mask].std(axis=1, ddof=1, keepdims=True)
        climatology[:, mask] = matrix[:, mask] - means
        standardised[:, mask] = climatology[:, mask] / standard_deviation
    for year in np.unique(years):
        for month in (6, 7, 8):
            mask = (years == year) & (months == month)
            within_record[:, mask] = matrix[:, mask] - matrix[:, mask].mean(
                axis=1, keepdims=True)
    return {
        "raw_wbt": matrix,
        "site_month_climatology_anomaly": climatology,
        "site_month_standardised_anomaly": standardised,
        "site_year_month_anomaly": within_record,
    }


def seasonal_progression_diagnostics(
        matrix: np.ndarray, day_rows: pd.DataFrame,
        operators: list[dict], labels: list[str]
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Check within-month timing and two calendar-progression adjustments."""
    timing_records = []
    for record_id, indices in day_rows.groupby(
            "record_id", sort=True).indices.items():
        subset = day_rows.iloc[indices]
        high = subset.regime.eq("high")
        middle = subset.regime.eq("middle")
        high_day = subset.loc[high, "date"].dt.day.mean()
        middle_day = subset.loc[middle, "date"].dt.day.mean()
        timing_records.append({
            "record_id": int(record_id),
            "year": int(record_id) // 100,
            "month": int(record_id) % 100,
            "high_mean_day_of_month": high_day,
            "middle_mean_day_of_month": middle_day,
            "high_minus_middle_day": high_day - middle_day,
        })
    timing_record = pd.DataFrame(timing_records)
    timing_year = timing_record.groupby("year", as_index=False).agg(
        high_mean_day_of_month=("high_mean_day_of_month", "mean"),
        middle_mean_day_of_month=("middle_mean_day_of_month", "mean"),
        high_minus_middle_day=("high_minus_middle_day", "mean"))
    timing_year["analysis_role"] = np.where(
        timing_year.year.isin(DISCOVERY), "development", "held_out")
    timing_summaries = []
    held = timing_record.loc[~timing_record.year.isin(DISCOVERY)]
    for month, group in [(0, held), *list(held.groupby("month", sort=True))]:
        if month == 0:
            yearly = group.groupby("year").high_minus_middle_day.mean()
            high_mean = group.high_mean_day_of_month.mean()
            middle_mean = group.middle_mean_day_of_month.mean()
            label = "June-August average"
        else:
            yearly = group.set_index("year").high_minus_middle_day
            high_mean = group.high_mean_day_of_month.mean()
            middle_mean = group.middle_mean_day_of_month.mean()
            label = str(int(month))
        timing_summaries.append({
            "month": label,
            "high_mean_day_of_month": high_mean,
            "middle_mean_day_of_month": middle_mean,
            **t_summary(yearly.to_numpy()),
        })
    timing_summary = pd.DataFrame(timing_summaries)

    # Remove a separate linear day-of-month trend at each site and record.
    linear_residual = np.zeros_like(matrix)
    for _, indices in day_rows.groupby("record_id", sort=True).indices.items():
        day = day_rows.iloc[indices].date.dt.day.to_numpy(dtype=float)
        centred_day = day - day.mean()
        fields = matrix[:, indices]
        beta = (fields @ centred_day) / np.square(centred_day).sum()
        fitted = fields.mean(axis=1, keepdims=True) + np.outer(
            beta, centred_day)
        linear_residual[:, indices] = fields - fitted

    # Subtract the site-specific calendar-day mean estimated without that year.
    month_day = day_rows.date.dt.strftime("%m-%d").to_numpy()
    years = day_rows.year.to_numpy()
    leave_one_year = np.zeros_like(matrix)
    for key in np.unique(month_day):
        key_indices = np.flatnonzero(month_day == key)
        fields = matrix[:, key_indices]
        if fields.shape[1] < 2:
            raise RuntimeError(f"Insufficient years for calendar day {key}")
        total = fields.sum(axis=1, keepdims=True)
        climatology_without_year = (total - fields) / (fields.shape[1] - 1)
        leave_one_year[:, key_indices] = fields - climatology_without_year

    records = pd.concat([
        record_contrasts(
            field_metrics(linear_residual, operators), day_rows, labels,
            "site_year_month_linear_detrended"),
        record_contrasts(
            field_metrics(leave_one_year, operators), day_rows, labels,
            "leave_one_year_daily_climatology_anomaly"),
    ], ignore_index=True)
    return timing_summary, timing_year, records


def hourly_sensitivity(labels: pd.DataFrame, operators: list[dict]) -> pd.DataFrame:
    records = []
    for path in sorted(HOURLY_DIR.glob("era5_land_*_jja_121sites.nc")):
        year = int(path.name.split("_")[2])
        dataset = xr.open_dataset(path, engine="h5netcdf")
        times = pd.DatetimeIndex(dataset.time.values)
        in_jja = times.month.isin([6, 7, 8])
        times = times[in_jja]
        temperature = dataset.t2m.values[:, in_jja]
        dewpoint = np.minimum(dataset.d2m.values[:, in_jja], temperature)
        pressure = dataset.sp.values[:, in_jja]
        wbt = BOLTON_WBT(temperature, dewpoint, pressure)
        metrics = field_metrics(wbt, operators)
        timestamp = pd.DataFrame({
            "time": times,
            "date": times.floor("D"),
            "hour": times.hour,
            "regional_mean_wbt": wbt.mean(axis=0),
        })
        timestamp = timestamp.merge(labels[["date", "record_id", "regime"]],
                                    on="date", how="left", validate="many_to_one")
        for hour in range(24):
            hour_indices = np.flatnonzero(timestamp.hour.to_numpy() == hour)
            hour_table = timestamp.iloc[hour_indices].copy().reset_index(drop=True)
            relabelled = classify_records(hour_table)
            for record_id, local_indices in hour_table.groupby(
                    "record_id", sort=True).indices.items():
                global_indices = hour_indices[np.asarray(local_indices)]
                for analysis, regime in (
                        ("fixed_hour", hour_table.iloc[local_indices].regime),
                        ("fixed_hour_relabelled",
                         relabelled.iloc[local_indices].regime)):
                    high = regime.to_numpy() == "high"
                    middle = regime.to_numpy() == "middle"
                    effects = metrics[global_indices][high].mean(axis=0) / (
                        metrics[global_indices][middle].mean(axis=0)) - 1
                    records.append({
                        "analysis": analysis,
                        "year": year,
                        "record_id": int(record_id),
                        "hour_utc": int(hour),
                        "profile_effect": effects.mean(),
                    })
        dates = timestamp.date.drop_duplicates().sort_values()
        date_codes = pd.Categorical(timestamp.date, categories=dates).codes
        daily_mean = np.column_stack([
            wbt[:, date_codes == code].mean(axis=1)
            for code in range(len(dates))
        ])
        daily_metrics = field_metrics(daily_mean, operators)
        daily_labels = labels.set_index("date").loc[dates].reset_index()
        daily_relabelled = daily_labels.copy()
        daily_relabelled["regional_mean_wbt"] = daily_mean.mean(axis=0)
        daily_relabelled = classify_records(daily_relabelled)
        for record_id, indices in daily_labels.groupby("record_id").indices.items():
            for analysis, regime in (
                    ("daily_mean_field", daily_labels.iloc[indices].regime),
                    ("daily_mean_relabelled",
                     daily_relabelled.iloc[indices].regime)):
                high = regime.to_numpy() == "high"
                middle = regime.to_numpy() == "middle"
                effects = daily_metrics[indices][high].mean(axis=0) / (
                    daily_metrics[indices][middle].mean(axis=0)) - 1
                records.append({
                    "analysis": analysis,
                    "year": year,
                    "record_id": int(record_id),
                    "hour_utc": -1,
                    "profile_effect": effects.mean(),
                })
        dataset.close()
        print(f"Hourly sensitivity: {year}")
    record = pd.DataFrame(records)
    year = record.groupby(["analysis", "year", "hour_utc"], as_index=False)[
        "profile_effect"].mean()
    held = year.loc[~year.year.isin(DISCOVERY)]
    summary = held.groupby(["analysis", "hour_utc"])["profile_effect"].apply(
        lambda x: pd.Series(t_summary(x.to_numpy()))).unstack().reset_index()
    return summary


def station_validation(
        event_manifest: pd.DataFrame, bandwidths: np.ndarray
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = PROJECT / "output_corrected" / "noaa_isd_era5_matched.csv.gz"
    data = pd.read_csv(path, parse_dates=["time_utc"])
    data["time_utc"] = pd.to_datetime(data.time_utc, utc=True).dt.tz_localize(None)
    data["date"] = data.time_utc.dt.floor("D")
    station = data.groupby("STATION").apply(
        lambda x: pd.Series({
            "n": len(x),
            "correlation": x.observed_wbt_c.corr(x.era5_wbt_c),
            "bias": (x.era5_wbt_c - x.observed_wbt_c).mean(),
            "mae": (x.era5_wbt_c - x.observed_wbt_c).abs().mean(),
            "rmse": np.sqrt(np.square(x.era5_wbt_c - x.observed_wbt_c).mean()),
        }), include_groups=False).reset_index()
    centred_observed = data.observed_wbt_c - data.groupby("STATION").observed_wbt_c.transform("mean")
    centred_era = data.era5_wbt_c - data.groupby("STATION").era5_wbt_c.transform("mean")
    summaries = pd.DataFrame([
        {
            "summary": "pooled",
            "n": len(data),
            "correlation": data.observed_wbt_c.corr(data.era5_wbt_c),
            "bias": (data.era5_wbt_c - data.observed_wbt_c).mean(),
            "mae": (data.era5_wbt_c - data.observed_wbt_c).abs().mean(),
            "rmse": np.sqrt(np.square(data.era5_wbt_c - data.observed_wbt_c).mean()),
        },
        {
            "summary": "within_station_centred",
            "n": len(data),
            "correlation": centred_observed.corr(centred_era),
        },
        {
            "summary": "station_equal_weight_mean",
            "n": len(station),
            "correlation": station.correlation.mean(),
            "bias": station.bias.mean(),
            "mae": station.mae.mean(),
            "rmse": station.rmse.mean(),
        },
        {
            "summary": "station_equal_weight_median",
            "n": len(station),
            "correlation": station.correlation.median(),
            "bias": station.bias.median(),
            "mae": station.mae.median(),
            "rmse": station.rmse.median(),
        },
    ])

    pair_parts = []
    keep = ["STATION", "time_utc", "LATITUDE", "LONGITUDE",
            "observed_wbt_c", "era5_wbt_c"]
    for _, group in data[keep].groupby("time_utc", sort=False):
        if len(group) < 2:
            continue
        group = group.sort_values("STATION")
        left, right = np.triu_indices(len(group), 1)
        values = group.reset_index(drop=True)
        lat1 = np.deg2rad(values.LATITUDE.to_numpy()[left])
        lat2 = np.deg2rad(values.LATITUDE.to_numpy()[right])
        dlat = lat2 - lat1
        dlon = np.deg2rad(values.LONGITUDE.to_numpy()[right] -
                          values.LONGITUDE.to_numpy()[left])
        haversine = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        distance = 6371.0088 * 2 * np.arcsin(np.sqrt(haversine))
        observed_difference = (values.observed_wbt_c.to_numpy()[left] -
                               values.observed_wbt_c.to_numpy()[right])
        era_difference = (values.era5_wbt_c.to_numpy()[left] -
                          values.era5_wbt_c.to_numpy()[right])
        pair_parts.append(pd.DataFrame({
            "time_utc": values.time_utc.iloc[0],
            "date": values.time_utc.iloc[0].floor("D"),
            "station_left": values.STATION.to_numpy()[left],
            "station_right": values.STATION.to_numpy()[right],
            "distance_km": distance,
            "observed_difference": observed_difference,
            "era_difference": era_difference,
            "observed_semivariance": np.square(observed_difference) / 2,
            "era_semivariance": np.square(era_difference) / 2,
        }))
    pairs = pd.concat(pair_parts, ignore_index=True)
    bins = [0, 250, 500, 1000, 2000, np.inf]
    labels = ["0-250", "250-500", "500-1000", "1000-2000", "2000+"]
    pairs["distance_bin_km"] = pd.cut(pairs.distance_km, bins=bins,
                                      labels=labels, right=False)
    pairs["observed_semivariance_residual"] = (
        pairs.observed_semivariance - pairs.groupby(
            "distance_bin_km", observed=True
        ).observed_semivariance.transform("mean"))
    pairs["era_semivariance_residual"] = (
        pairs.era_semivariance - pairs.groupby(
            "distance_bin_km", observed=True
        ).era_semivariance.transform("mean"))
    spatial = pairs.groupby("distance_bin_km", observed=True).agg(
        pairs=("distance_km", "size"),
        distance_mean_km=("distance_km", "mean"),
        observed_semivariance=("observed_semivariance", "mean"),
        era_semivariance=("era_semivariance", "mean"),
    ).reset_index()
    within_bin = pairs.groupby("distance_bin_km", observed=True).apply(
        lambda x: pd.Series({
            "within_bin_difference_correlation":
                x.observed_difference.corr(x.era_difference),
            "within_bin_semivariance_correlation":
                x.observed_semivariance.corr(x.era_semivariance),
        }), include_groups=False).reset_index()
    spatial = spatial.merge(within_bin, on="distance_bin_km",
                            validate="one_to_one")
    spatial["signed_pair_difference_correlation"] = pairs.observed_difference.corr(
        pairs.era_difference)
    spatial["pair_semivariance_correlation"] = pairs.observed_semivariance.corr(
        pairs.era_semivariance)
    spatial["distance_bin_residual_semivariance_correlation"] = (
        pairs.observed_semivariance_residual.corr(
            pairs.era_semivariance_residual))

    equal_pairs = pairs.groupby(
        ["station_left", "station_right"], as_index=False).agg(
            hours=("time_utc", "size"),
            distance_km=("distance_km", "mean"),
            observed_semivariance=("observed_semivariance", "mean"),
            era_semivariance=("era_semivariance", "mean"))
    equal_pairs["distance_bin_km"] = pd.cut(
        equal_pairs.distance_km, bins=bins, labels=labels, right=False)
    equal_pairs["observed_residual"] = (
        equal_pairs.observed_semivariance - equal_pairs.groupby(
            "distance_bin_km", observed=True
        ).observed_semivariance.transform("mean"))
    equal_pairs["era_residual"] = (
        equal_pairs.era_semivariance - equal_pairs.groupby(
            "distance_bin_km", observed=True
        ).era_semivariance.transform("mean"))

    pair_temporal = pairs.groupby(
        ["station_left", "station_right"]).apply(
            lambda x: pd.Series({
                "hours": len(x),
                "difference_correlation":
                    x.observed_difference.corr(x.era_difference),
                "semivariance_correlation":
                    x.observed_semivariance.corr(x.era_semivariance),
            }), include_groups=False).reset_index()

    # Date-block bootstrap for the distance-bin-centred semivariance correlation.
    block = pairs.groupby("date").apply(
        lambda x: pd.Series({
            "n": len(x),
            "sum_x": x.observed_semivariance_residual.sum(),
            "sum_y": x.era_semivariance_residual.sum(),
            "sum_x2": np.square(x.observed_semivariance_residual).sum(),
            "sum_y2": np.square(x.era_semivariance_residual).sum(),
            "sum_xy": (x.observed_semivariance_residual *
                       x.era_semivariance_residual).sum(),
        }), include_groups=False).reset_index()
    sufficient = block[["n", "sum_x", "sum_y", "sum_x2", "sum_y2",
                        "sum_xy"]].to_numpy()

    def correlation_from_sufficient(values: np.ndarray) -> float:
        n, sum_x, sum_y, sum_x2, sum_y2, sum_xy = values.sum(axis=0)
        covariance = sum_xy - sum_x * sum_y / n
        variance_x = sum_x2 - sum_x**2 / n
        variance_y = sum_y2 - sum_y**2 / n
        return float(covariance / math.sqrt(variance_x * variance_y))

    rng = np.random.default_rng(20260806)
    bootstrap = np.empty(1999)
    for b in range(len(bootstrap)):
        sampled = rng.integers(0, len(sufficient), size=len(sufficient))
        bootstrap[b] = correlation_from_sufficient(sufficient[sampled])

    advanced = pd.DataFrame([
        {
            "diagnostic": "pair_hour_semivariance_correlation",
            "estimate": pairs.observed_semivariance.corr(
                pairs.era_semivariance),
            "n": len(pairs),
        },
        {
            "diagnostic": "distance_bin_residual_semivariance_correlation",
            "estimate": pairs.observed_semivariance_residual.corr(
                pairs.era_semivariance_residual),
            "ci_lower": np.quantile(bootstrap, 0.025),
            "ci_upper": np.quantile(bootstrap, 0.975),
            "n": len(pairs),
            "clusters": len(block),
        },
        {
            "diagnostic": "station_pair_equal_weight_semivariance_correlation",
            "estimate": equal_pairs.observed_semivariance.corr(
                equal_pairs.era_semivariance),
            "n": len(equal_pairs),
        },
        {
            "diagnostic":
                "station_pair_equal_weight_residual_semivariance_correlation",
            "estimate": equal_pairs.observed_residual.corr(
                equal_pairs.era_residual),
            "n": len(equal_pairs),
        },
        {
            "diagnostic": "mean_within_station_pair_difference_correlation",
            "estimate": pair_temporal.difference_correlation.mean(),
            "ci_lower": pair_temporal.difference_correlation.quantile(0.025),
            "ci_upper": pair_temporal.difference_correlation.quantile(0.975),
            "n": pair_temporal.difference_correlation.notna().sum(),
        },
        {
            "diagnostic": "mean_within_station_pair_semivariance_correlation",
            "estimate": pair_temporal.semivariance_correlation.mean(),
            "ci_lower": pair_temporal.semivariance_correlation.quantile(0.025),
            "ci_upper": pair_temporal.semivariance_correlation.quantile(0.975),
            "n": pair_temporal.semivariance_correlation.notna().sum(),
        },
    ])

    # Compare graph dispersion on the same station subset at every matched hour.
    graph_rows = []
    for time, group in data.groupby("time_utc", sort=False):
        group = group.sort_values("STATION").drop_duplicates("STATION")
        if len(group) < 10:
            continue
        latitude = np.deg2rad(group.LATITUDE.to_numpy())
        longitude = np.deg2rad(group.LONGITUDE.to_numpy())
        dlat = latitude[:, None] - latitude[None, :]
        dlon = longitude[:, None] - longitude[None, :]
        haversine = (np.sin(dlat / 2)**2 +
                     np.cos(latitude[:, None]) *
                     np.cos(latitude[None, :]) * np.sin(dlon / 2)**2)
        distance_matrix = 6371.0088 * 2 * np.arcsin(
            np.sqrt(np.clip(haversine, 0, 1)))
        upper = np.triu_indices(len(group), 1)
        observed = group.observed_wbt_c.to_numpy()
        era = group.era5_wbt_c.to_numpy()
        for bandwidth in bandwidths:
            weight = np.exp(-np.square(distance_matrix[upper]) /
                            (2 * bandwidth**2))
            observed_q = np.sum(
                weight * np.square(observed[upper[0]] - observed[upper[1]])) / (
                    2 * weight.sum())
            era_q = np.sum(
                weight * np.square(era[upper[0]] - era[upper[1]])) / (
                    2 * weight.sum())
            graph_rows.append({
                "time_utc": time,
                "bandwidth_km": bandwidth,
                "stations": len(group),
                "observed_q": observed_q,
                "era_q": era_q,
            })
    graph = pd.DataFrame(graph_rows)
    event = event_manifest.loc[
        event_manifest.year.isin(DISCOVERY),
        ["peak_time", "record_id", "year", "month", "regime"]].copy()
    event["peak_time"] = pd.to_datetime(event.peak_time).dt.tz_localize(None)
    graph_event = graph.merge(
        event, left_on="time_utc", right_on="peak_time", how="inner",
        validate="many_to_one")
    contrast_rows = []
    for (record_id, bandwidth), group in graph_event.groupby(
            ["record_id", "bandwidth_km"], sort=True):
        high = group.regime.eq("high")
        middle = group.regime.eq("middle")
        if high.sum() < 2 or middle.sum() < 2:
            continue
        contrast_rows.append({
            "record_id": int(record_id),
            "year": int(record_id) // 100,
            "month": int(record_id) % 100,
            "bandwidth_km": bandwidth,
            "high_fields": int(high.sum()),
            "middle_fields": int(middle.sum()),
            "observed_effect": group.loc[high, "observed_q"].mean() /
                               group.loc[middle, "observed_q"].mean() - 1,
            "era_effect": group.loc[high, "era_q"].mean() /
                          group.loc[middle, "era_q"].mean() - 1,
        })
    contrasts = pd.DataFrame(contrast_rows)
    graph_agreement = graph.groupby("bandwidth_km").apply(
        lambda x: pd.Series({
            "analysis": "matched_hour_graph_dispersion",
            "fields": len(x),
            "records": np.nan,
            "observed_value": x.observed_q.mean(),
            "era_value": x.era_q.mean(),
            "correlation": x.observed_q.corr(x.era_q),
        }), include_groups=False).reset_index()
    event_profile = contrasts.groupby("bandwidth_km").apply(
        lambda x: pd.Series({
            "analysis": "development_event_contrast",
            "fields": int(x.high_fields.sum() + x.middle_fields.sum()),
            "records": len(x),
            "observed_value": x.observed_effect.mean(),
            "era_value": x.era_effect.mean(),
            "correlation": np.nan,
        }), include_groups=False).reset_index()
    graph_validation = pd.concat(
        [graph_agreement, event_profile], ignore_index=True)
    return summaries, spatial, advanced, graph_validation


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(DAILY_DIR.glob("era5_land_*_jja_daily_fields.csv.gz"))
    columns = ["year", "month", "record_id", "analysis_date", "time",
               "day_definition", "site_id", "requested_lon", "requested_lat",
               "regional_mean_wbt", "wbt"]
    daily = pd.concat([pd.read_csv(path, usecols=columns) for path in paths],
                      ignore_index=True)
    daily = daily.loc[daily.day_definition.eq("utc")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    daily["time"] = pd.to_datetime(daily.time)
    sites = daily[["site_id", "requested_lon", "requested_lat"]].drop_duplicates(
    ).sort_values("site_id").reset_index(drop=True)
    distance = projected_distance(sites)
    gaussian = gaussian_operators(distance)

    day_rows = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"),
        peak_time=("time", "first"))
    day_rows = classify_records(day_rows).sort_values(
        ["record_id", "date"]).reset_index(drop=True)
    day_rows.attrs["sites"] = sites
    event_manifest = day_rows.copy()
    event_manifest["quantile_type"] = "linear interpolation (R type 7)"
    event_manifest["ties_at_q25"] = np.isclose(event_manifest.regional_mean_wbt,
                                                event_manifest.q25)
    event_manifest["ties_at_q75"] = np.isclose(event_manifest.regional_mean_wbt,
                                                event_manifest.q75)
    event_manifest.to_csv(OUTPUT / "sensitivity_event_manifest.csv", index=False)

    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    pivot = pivot.reindex(index=sites.site_id, columns=day_rows.date)
    matrix = pivot.to_numpy()
    if not np.isfinite(matrix).all():
        raise RuntimeError("Incomplete daily WBT matrix")

    labels = [f"{operator['bandwidth_km']:.6f}" for operator in gaussian]
    raw_metrics = field_metrics(matrix, gaussian)
    raw_records = record_contrasts(raw_metrics, day_rows, labels, "equal_site_ratio")

    # Area-weighted graph dispersion, retaining the prespecified event labels.
    area = np.cos(np.deg2rad(sites.requested_lat.to_numpy()))
    area_operators = []
    for base in gaussian:
        weighted = base["W"] * np.outer(area, area)
        operator = graph_operator(weighted)
        operator.update(h_factor=base["h_factor"], bandwidth_km=base["bandwidth_km"],
                        kernel="Gaussian area-weighted")
        area_operators.append(operator)
    area_records = record_contrasts(field_metrics(matrix, area_operators), day_rows,
                                    labels, "area_weighted_fixed_labels")
    area_day_rows = day_rows.copy()
    area_day_rows["regional_mean_wbt"] = np.average(
        matrix, axis=0, weights=area)
    area_day_rows = classify_records(area_day_rows)
    area_records_relabelled = record_contrasts(
        field_metrics(matrix, area_operators), area_day_rows, labels,
        "area_weighted_relabelled_primary_hours")

    # Kernel sensitivity matched by weighted mean pair distance.
    kernel_records = [raw_records, area_records, area_records_relabelled]
    matched_operators = []
    gaussian_metadata = operator_metadata(gaussian, distance, len(sites))
    for kernel in ("Exponential", "Compact quadratic"):
        operators = []
        for j, base in enumerate(gaussian):
            weights, parameter = matched_kernel(
                distance, gaussian_metadata.loc[j, "weighted_distance_mean_km"], kernel)
            operator = graph_operator(weights)
            operator.update(h_factor=base["h_factor"], bandwidth_km=parameter,
                            kernel=kernel)
            operators.append(operator)
        matched_operators.extend(operators)
        kernel_records.append(record_contrasts(
            field_metrics(matrix, operators), day_rows, labels,
            f"{kernel.lower().replace(' ', '_')}_matched_distance"))
    all_records = pd.concat(kernel_records, ignore_index=True)

    # Raw-field and anomaly-field sensitivity.
    for name, transformed in anomaly_matrices(matrix, day_rows, sites).items():
        if name == "raw_wbt":
            continue
        all_records = pd.concat([
            all_records,
            record_contrasts(field_metrics(transformed, gaussian), day_rows,
                             labels, name),
        ], ignore_index=True)

    timing_summary, timing_year, seasonal_records = (
        seasonal_progression_diagnostics(
            matrix, day_rows, gaussian, labels))
    all_records = pd.concat(
        [all_records, seasonal_records], ignore_index=True)
    timing_summary.to_csv(
        OUTPUT / "sensitivity_seasonal_timing.csv", index=False)
    timing_year.to_csv(
        OUTPUT / "sensitivity_seasonal_timing_yearly.csv", index=False)

    robustness = pd.concat([
        summarize_analysis(all_records, "ratio_effect"),
        summarize_analysis(raw_records, "log_effect"),
    ], ignore_index=True)
    robustness.to_csv(OUTPUT / "sensitivity_robustness_summary.csv", index=False)
    summarize_scales(all_records).to_csv(
        OUTPUT / "sensitivity_scale_sensitivity.csv", index=False)

    denominator = raw_records.loc[~raw_records.year.isin(DISCOVERY)].groupby(
        "scale").mean_middle.agg([
            "count", "min", lambda x: x.quantile(.01), lambda x: x.quantile(.05),
            lambda x: x.quantile(.25), "median", lambda x: x.quantile(.75),
            "max", "mean", "std"])
    denominator.columns = ["records", "min", "q01", "q05", "q25", "median",
                           "q75", "max", "mean", "sd"]
    denominator["coefficient_of_variation"] = denominator.sd / denominator["mean"]
    denominator = denominator.reset_index()
    metadata = pd.concat([
        gaussian_metadata,
        operator_metadata(matched_operators, distance, len(sites)),
    ], ignore_index=True)
    denominator.merge(gaussian_metadata, left_on="scale",
                      right_on=gaussian_metadata.bandwidth_km.map(
                          lambda x: f"{x:.6f}"), how="left").to_csv(
        OUTPUT / "sensitivity_scale_diagnostics.csv", index=False)

    raw_yearly = year_effects(raw_records, "ratio_effect").rename(
        columns={"ratio_effect": "ratio_effect"})
    temporal_diagnostics(raw_yearly).to_csv(
        OUTPUT / "sensitivity_temporal_diagnostics.csv", index=False)
    gradient_summary, gradient_yearly = gradient_decomposition(
        matrix, day_rows, gaussian)
    gradient_summary.to_csv(
        OUTPUT / "sensitivity_gradient_decomposition.csv", index=False)
    gradient_yearly.to_csv(
        OUTPUT / "sensitivity_gradient_yearly.csv", index=False)
    metadata.to_csv(OUTPUT / "sensitivity_kernel_metadata.csv", index=False)

    hour = hourly_sensitivity(day_rows[["date", "record_id", "regime"]], gaussian)
    hour.to_csv(OUTPUT / "sensitivity_hourly_sensitivity.csv", index=False)

    marginal, spatial, station_advanced, station_graph = station_validation(
        event_manifest, np.array(
            [operator["bandwidth_km"] for operator in gaussian]))
    marginal.to_csv(OUTPUT / "sensitivity_station_validation.csv", index=False)
    spatial.to_csv(OUTPUT / "sensitivity_station_spatial_validation.csv", index=False)
    station_advanced.to_csv(
        OUTPUT / "sensitivity_station_advanced_validation.csv", index=False)
    station_graph.to_csv(
        OUTPUT / "sensitivity_station_graph_validation.csv", index=False)

    print("sensitivity diagnostics complete")
    print(robustness.to_string(index=False))
    print(temporal_diagnostics(raw_yearly).to_string(index=False))


if __name__ == "__main__":
    main()
