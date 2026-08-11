#!/usr/bin/env python3
"""Post-analysis continuous-intensity and calendar-matching diagnostics.

This script implements the definitions frozen in EXTENSION_ANALYSIS_PROTOCOL.md.
It does not alter the primary quartile estimator or any manuscript source.
Only the original UTC daily fields and the 121-site graph are used.  The two
development summers (2015 and 2022) are retained in record-level output but
excluded from the reported 33-summer summaries.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
DAILY_DIR = PROJECT / "data" / "era5_confirmatory" / "daily_fields"
PRIMARY_OUTPUT = PROJECT / "output_confirmatory"
OUTPUT = PROJECT / "output_extension_methods"
DISCOVERY_YEARS = {2015, 2022}
WINDOWS = (3, 5)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def projected_distance(sites: pd.DataFrame) -> np.ndarray:
    lat0 = np.deg2rad(sites.requested_lat.mean())
    coordinates = np.column_stack((
        sites.requested_lon.to_numpy() * 111.32 * np.cos(lat0),
        sites.requested_lat.to_numpy() * 110.57,
    ))
    differences = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.square(differences).sum(axis=2))


def edge_metrics(matrix: np.ndarray, distance: np.ndarray,
                 bandwidths: np.ndarray) -> np.ndarray:
    upper = np.triu_indices_from(distance, 1)
    weights = np.exp(
        -np.square(distance[upper][:, None]) /
        (2 * np.square(bandwidths[None, :]))
    )
    weight_sum = weights.sum(axis=0)
    result = np.empty((matrix.shape[1], len(bandwidths)), dtype=float)
    for start in range(0, matrix.shape[1], 128):
        stop = min(start + 128, matrix.shape[1])
        differences = (matrix[upper[0], start:stop] -
                       matrix[upper[1], start:stop])
        result[start:stop] = np.square(differences).T @ weights / (
            2 * weight_sum)
    return result


def load_fields() -> tuple[pd.DataFrame, np.ndarray, np.ndarray, list[Path]]:
    paths = sorted(DAILY_DIR.glob("era5_land_*_jja_daily_fields.csv.gz"))
    if len(paths) != 35:
        raise RuntimeError(f"Expected 35 daily-field files, found {len(paths)}")
    columns = [
        "year", "month", "record_id", "analysis_date", "day_definition",
        "site_id", "requested_lon", "requested_lat", "regional_mean_wbt",
        "wbt",
    ]
    daily = pd.concat(
        [pd.read_csv(path, usecols=columns) for path in paths],
        ignore_index=True,
    )
    daily = daily.loc[daily.day_definition.eq("utc")].copy()
    daily["date"] = pd.to_datetime(daily.analysis_date)
    sites = (daily[["site_id", "requested_lon", "requested_lat"]]
             .drop_duplicates().sort_values("site_id").reset_index(drop=True))
    if len(sites) != 121:
        raise RuntimeError(f"Expected 121 sites, found {len(sites)}")

    day_rows = daily.groupby(["record_id", "date"], as_index=False).agg(
        year=("year", "first"), month=("month", "first"),
        regional_mean_wbt=("regional_mean_wbt", "first"),
    ).sort_values(["record_id", "date"]).reset_index(drop=True)
    if len(day_rows) != 35 * 92:
        raise RuntimeError(f"Expected 3,220 UTC days, found {len(day_rows)}")

    quantiles = day_rows.groupby("record_id").regional_mean_wbt.quantile(
        [0.25, 0.75], interpolation="linear").unstack()
    quantiles.columns = ["q25", "q75"]
    day_rows = day_rows.merge(
        quantiles, left_on="record_id", right_index=True,
        validate="many_to_one",
    )
    day_rows["regime"] = np.where(
        day_rows.regional_mean_wbt >= day_rows.q75, "high",
        np.where(day_rows.regional_mean_wbt >= day_rows.q25, "middle", "low"),
    )

    pivot = daily.pivot(index="site_id", columns="date", values="wbt")
    pivot = pivot.reindex(index=sites.site_id, columns=day_rows.date)
    matrix = pivot.to_numpy()
    if matrix.shape != (121, 35 * 92) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid site-by-day matrix {matrix.shape}")

    metadata = pd.read_csv(PRIMARY_OUTPUT / "confirmatory_graph_metadata.csv")
    bandwidths = (metadata.loc[metadata.definition_index.eq(1), "bandwidth_km"]
                  .to_numpy(dtype=float))
    if len(bandwidths) != 5:
        raise RuntimeError("Expected five fixed bandwidths")
    metrics = edge_metrics(matrix, projected_distance(sites), bandwidths)
    if not np.isfinite(metrics).all() or (metrics <= 0).any():
        raise RuntimeError("Graph dispersions must be finite and positive")
    return day_rows, metrics, bandwidths, paths


def summarize(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    estimate = float(values.mean())
    se = float(values.std(ddof=1) / math.sqrt(n))
    critical = float(stats.t.ppf(0.975, n - 1))
    return {
        "summers": n,
        "estimate": estimate,
        "standard_error": se,
        "ci_lower": estimate - critical * se,
        "ci_upper": estimate + critical * se,
        "negative_summers": int((values < 0).sum()),
    }


def continuous_intensity(day_rows: pd.DataFrame, metrics: np.ndarray,
                         bandwidths: np.ndarray) -> tuple[pd.DataFrame, ...]:
    monthly_rows: list[dict] = []
    for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
        indices = np.asarray(indices, dtype=int)
        subset = day_rows.iloc[indices]
        x = subset.regional_mean_wbt.to_numpy(dtype=float)
        x -= x.mean()
        denominator = float(x @ x)
        if denominator <= 0:
            raise RuntimeError(f"Zero WBT variance in record {record_id}")
        log_q = np.log(metrics[indices])
        slopes = x @ log_q / denominator
        year = int(subset.year.iloc[0])
        month = int(subset.month.iloc[0])
        for bandwidth, slope in zip(bandwidths, slopes):
            monthly_rows.append({
                "record_id": int(record_id), "year": year, "month": month,
                "analysis_role": "discovery" if year in DISCOVERY_YEARS
                else "held_out", "bandwidth_km": bandwidth,
                "slope_log_q_per_c": float(slope), "days": len(indices),
            })
    monthly = pd.DataFrame(monthly_rows)
    yearly = monthly.groupby(
        ["year", "analysis_role", "bandwidth_km"], as_index=False,
    ).agg(slope_log_q_per_c=("slope_log_q_per_c", "mean"),
          months=("month", "nunique"))
    if not yearly.months.eq(3).all():
        raise RuntimeError("Continuous slope hierarchy has incomplete summers")
    held = yearly.loc[yearly.analysis_role.eq("held_out")]
    scale_rows = []
    for bandwidth, group in held.groupby("bandwidth_km", sort=True):
        scale_rows.append({"bandwidth_km": bandwidth,
                           **summarize(group.slope_log_q_per_c.to_numpy())})
    scale_summary = pd.DataFrame(scale_rows)
    profile_year = held.groupby("year", as_index=False).agg(
        slope_log_q_per_c=("slope_log_q_per_c", "mean"),
        bandwidths=("bandwidth_km", "nunique"),
    )
    if not profile_year.bandwidths.eq(5).all():
        raise RuntimeError("Continuous profile has incomplete bandwidths")
    profile_summary = pd.DataFrame([summarize(
        profile_year.slope_log_q_per_c.to_numpy())])
    return monthly, yearly, scale_summary, profile_year, profile_summary


def calendar_matching(day_rows: pd.DataFrame, metrics: np.ndarray,
                      bandwidths: np.ndarray) -> tuple[pd.DataFrame, ...]:
    record_rows: list[dict] = []
    for window in WINDOWS:
        for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
            indices = np.asarray(indices, dtype=int)
            subset = day_rows.iloc[indices]
            high_local = np.flatnonzero(subset.regime.to_numpy() == "high")
            middle_local = np.flatnonzero(subset.regime.to_numpy() == "middle")
            dates = subset.date.to_numpy(dtype="datetime64[D]")
            high_effects = []
            matched_middle_counts = []
            for high_index in high_local:
                differences = np.abs(
                    (dates[middle_local] - dates[high_index]).astype(int))
                eligible = middle_local[differences <= window]
                if not len(eligible):
                    continue
                denominator = metrics[indices[eligible]].mean(axis=0)
                high_effects.append(metrics[indices[high_index]] / denominator - 1)
                matched_middle_counts.append(len(eligible))
            year = int(subset.year.iloc[0])
            month = int(subset.month.iloc[0])
            available = len(high_effects)
            effects = (np.vstack(high_effects).mean(axis=0) if available else
                       np.repeat(np.nan, len(bandwidths)))
            for bandwidth, effect in zip(bandwidths, effects):
                record_rows.append({
                    "window_days": window, "record_id": int(record_id),
                    "year": year, "month": month,
                    "analysis_role": "discovery" if year in DISCOVERY_YEARS
                    else "held_out", "bandwidth_km": bandwidth,
                    "matched_effect": effect,
                    "high_days_total": len(high_local),
                    "high_days_matched": available,
                    "high_day_coverage": available / len(high_local),
                    "mean_eligible_middle_days": (float(np.mean(matched_middle_counts))
                                                  if available else 0.0),
                    "record_has_no_match": available == 0,
                })
    records = pd.DataFrame(record_rows)
    if records.matched_effect.isna().any():
        # Records with no eligible match are reported, but cannot enter the
        # equal-record estimator.  No such records are expected for these data.
        bad = records.loc[records.matched_effect.isna(),
                          ["window_days", "record_id"]].drop_duplicates()
        raise RuntimeError(f"Calendar match has empty records:\n{bad}")
    yearly = records.groupby(
        ["window_days", "year", "analysis_role", "bandwidth_km"],
        as_index=False,
    ).agg(matched_effect=("matched_effect", "mean"),
          months=("month", "nunique"))
    if not yearly.months.eq(3).all():
        raise RuntimeError("Calendar-matching hierarchy has incomplete summers")
    held = yearly.loc[yearly.analysis_role.eq("held_out")]
    scale_rows = []
    for (window, bandwidth), group in held.groupby(
            ["window_days", "bandwidth_km"], sort=True):
        scale_rows.append({"window_days": window, "bandwidth_km": bandwidth,
                           **summarize(group.matched_effect.to_numpy())})
    scale_summary = pd.DataFrame(scale_rows)
    profile_year = held.groupby(["window_days", "year"], as_index=False).agg(
        matched_effect=("matched_effect", "mean"),
        bandwidths=("bandwidth_km", "nunique"),
    )
    if not profile_year.bandwidths.eq(5).all():
        raise RuntimeError("Calendar profile has incomplete bandwidths")
    profile_rows = []
    for window, group in profile_year.groupby("window_days", sort=True):
        profile_rows.append({"window_days": window,
                             **summarize(group.matched_effect.to_numpy())})
    profile_summary = pd.DataFrame(profile_rows)
    coverage = (records.loc[records.analysis_role.eq("held_out")]
                .drop_duplicates(["window_days", "record_id"]))
    coverage_summary = coverage.groupby("window_days", as_index=False).agg(
        records=("record_id", "nunique"),
        high_days_total=("high_days_total", "sum"),
        high_days_matched=("high_days_matched", "sum"),
        records_with_no_match=("record_has_no_match", "sum"),
        mean_record_coverage=("high_day_coverage", "mean"),
        mean_eligible_middle_days=("mean_eligible_middle_days", "mean"),
    )
    coverage_summary["pooled_high_day_coverage"] = (
        coverage_summary.high_days_matched / coverage_summary.high_days_total)
    return (records, yearly, scale_summary, profile_year, profile_summary,
            coverage_summary)


def write_tex(continuous_scale: pd.DataFrame,
              continuous_profile: pd.DataFrame,
              matching_scale: pd.DataFrame,
              matching_profile: pd.DataFrame,
              coverage: pd.DataFrame) -> list[Path]:
    continuous_path = OUTPUT / "extension_continuous_intensity_table.tex"
    lines = [
        "% Generated by code/44_extension_empirical_methods.py",
        "\\begin{tabular}{rrrrr}",
        "\\toprule",
        "Bandwidth (km) & Slope ($^\\circ$C$^{-1}$) & 95\\% CI & Summers & Negative \\\\",
        "\\midrule",
    ]
    for row in continuous_scale.itertuples(index=False):
        lines.append(
            f"{row.bandwidth_km:,.0f} & {row.estimate:.4f} & "
            f"[{row.ci_lower:.4f}, {row.ci_upper:.4f}] & "
            f"{row.summers} & {row.negative_summers} \\\\"
        )
    p = continuous_profile.iloc[0]
    lines += ["\\midrule",
              f"Five-scale mean & {p.estimate:.4f} & "
              f"[{p.ci_lower:.4f}, {p.ci_upper:.4f}] & "
              f"{int(p.summers)} & {int(p.negative_summers)} \\\\ ",
              "\\bottomrule", "\\end{tabular}", ""]
    continuous_path.write_text("\n".join(lines), encoding="utf-8")

    matching_path = OUTPUT / "extension_calendar_matching_table.tex"
    lines = [
        "% Generated by code/44_extension_empirical_methods.py",
        "\\begin{tabular}{rrrrrr}",
        "\\toprule",
        "Window & Bandwidth (km) & Effect (\\%) & 95\\% CI (\\%) & Summers & Negative \\\\",
        "\\midrule",
    ]
    for window in WINDOWS:
        group = matching_scale.loc[matching_scale.window_days.eq(window)]
        for row in group.itertuples(index=False):
            lines.append(
                f"$\\pm${window} d & {row.bandwidth_km:,.0f} & "
                f"{100 * row.estimate:.2f} & "
                f"[{100 * row.ci_lower:.2f}, {100 * row.ci_upper:.2f}] & "
                f"{row.summers} & {row.negative_summers} \\\\"
            )
        p = matching_profile.loc[matching_profile.window_days.eq(window)].iloc[0]
        c = coverage.loc[coverage.window_days.eq(window)].iloc[0]
        lines.append(
            f"$\\pm${window} d & Five-scale mean & {100 * p.estimate:.2f} & "
            f"[{100 * p.ci_lower:.2f}, {100 * p.ci_upper:.2f}] & "
            f"{int(p.summers)} & {int(p.negative_summers)} \\\\"
        )
        lines.append(
            f"\\multicolumn{{6}}{{l}}{{\\footnotesize Coverage: "
            f"{int(c.high_days_matched)}/{int(c.high_days_total)} high days "
            f"({100 * c.pooled_high_day_coverage:.1f}\\%); "
            f"{int(c.records_with_no_match)} records without a match.}} \\\\"
        )
        if window != WINDOWS[-1]:
            lines.append("\\addlinespace")
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    matching_path.write_text("\n".join(lines), encoding="utf-8")
    return [continuous_path, matching_path]


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    day_rows, metrics, bandwidths, inputs = load_fields()

    # Reproduce the frozen raw-ratio profile before running new diagnostics.
    raw_records = []
    for record_id, indices in day_rows.groupby("record_id", sort=True).indices.items():
        indices = np.asarray(indices, dtype=int)
        subset = day_rows.iloc[indices]
        high = subset.regime.to_numpy() == "high"
        middle = subset.regime.to_numpy() == "middle"
        ratio = (metrics[indices][high].mean(axis=0) /
                 metrics[indices][middle].mean(axis=0) - 1)
        for bandwidth, effect in zip(bandwidths, ratio):
            raw_records.append({"year": int(subset.year.iloc[0]),
                                "bandwidth_km": bandwidth, "effect": effect})
    check = pd.DataFrame(raw_records)
    reproduced = (check.loc[~check.year.isin(DISCOVERY_YEARS)]
                  .groupby(["year", "bandwidth_km"], as_index=False).effect.mean()
                  .groupby("year").effect.mean().mean())
    expected = pd.read_csv(PRIMARY_OUTPUT / "confirmatory_primary_results.csv")
    expected = float(expected.loc[expected.day_definition.eq("utc"), "estimate"].iloc[0])
    if abs(reproduced - expected) > 1e-10:
        raise RuntimeError(f"Primary reproduction failed: {reproduced} vs {expected}")

    continuous = continuous_intensity(day_rows, metrics, bandwidths)
    matching = calendar_matching(day_rows, metrics, bandwidths)
    continuous_names = [
        "extension_continuous_month_slopes.csv",
        "extension_continuous_year_scale_slopes.csv",
        "extension_continuous_scale_summary.csv",
        "extension_continuous_year_profile_slopes.csv",
        "extension_continuous_profile_summary.csv",
    ]
    matching_names = [
        "extension_calendar_matching_records.csv",
        "extension_calendar_matching_year_scale.csv",
        "extension_calendar_matching_scale_summary.csv",
        "extension_calendar_matching_year_profile.csv",
        "extension_calendar_matching_profile_summary.csv",
        "extension_calendar_matching_coverage.csv",
    ]
    output_paths = []
    for table, name in zip(continuous, continuous_names):
        path = OUTPUT / name
        table.to_csv(path, index=False)
        output_paths.append(path)
    for table, name in zip(matching, matching_names):
        path = OUTPUT / name
        table.to_csv(path, index=False)
        output_paths.append(path)
    tex_paths = write_tex(continuous[2], continuous[4], matching[2],
                          matching[4], matching[5])
    output_paths.extend(tex_paths)

    input_manifest = [{"path": str(path.relative_to(PROJECT)),
                       "bytes": path.stat().st_size, "sha256": sha256(path)}
                      for path in inputs]
    audit = {
        "protocol": "EXTENSION_ANALYSIS_PROTOCOL.md",
        "script": str(Path(__file__).resolve().relative_to(PROJECT)),
        "analysis_role": "post-analysis structural sensitivity",
        "random_seed": None,
        "day_definition": "UTC regional-mean peak",
        "site_count": 121,
        "summer_count_total": 35,
        "summer_count_held_out": 33,
        "month_year_records_total": 105,
        "fixed_bandwidths_km": bandwidths.tolist(),
        "primary_reproduction": {"computed": reproduced, "expected": expected,
                                 "absolute_error": abs(reproduced - expected)},
        "calendar_windows_days": list(WINDOWS),
        "input_files": input_manifest,
        "output_files": [
            {"path": str(path.relative_to(PROJECT)), "bytes": path.stat().st_size,
             "sha256": sha256(path)} for path in output_paths
        ],
    }
    audit_path = OUTPUT / "extension_empirical_methods_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"Primary reproduction: {reproduced:.12f}")
    print("Continuous five-scale summary:")
    print(continuous[4].to_string(index=False))
    print("Calendar-matching summaries:")
    print(matching[4].to_string(index=False))
    print("Coverage:")
    print(matching[5].to_string(index=False))


if __name__ == "__main__":
    main()
