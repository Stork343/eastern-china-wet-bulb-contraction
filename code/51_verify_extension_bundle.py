#!/usr/bin/env python3
"""Strictly verify every artifact in the extension analysis bundle.

This verifier is intentionally fail-closed.  In particular, an unfinished NOAA
download, ERA5-Land station extraction, or NOAA--ERA5-Land analysis is an error;
no stage is treated as optional.  The acquisition audit also preserves the
outcome-blind administrative operationalisation: the frozen END>=2025-08-31
filter returned zero stations because the official in-scope history snapshot
ended on 2025-08-24, so only that administrative boundary was truncated.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


PROJECT = Path(__file__).resolve().parent.parent
HIST_DATA = PROJECT / "data" / "era5_historical_extension"
HIST_OUT = PROJECT / "output_historical_extension"
METHOD_OUT = PROJECT / "output_extension_methods"
ELEV_OUT = PROJECT / "output_elevation_basis"
NOAA_DATA = PROJECT / "data" / "noaa_isd_extension"
NOAA_ERA = NOAA_DATA / "era5_land_points"
NOAA_OUT = PROJECT / "output_noaa_extension"
PRIMARY_SITE_FILE = PROJECT / "data" / "grid" / "eastern_china_121_sites.csv"
MAIN_TEX = PROJECT / "manuscript" / "main.tex"
SUPPLEMENT_TEX = PROJECT / "manuscript" / "supplement_theory.tex"

HISTORICAL_YEARS = tuple(range(1950, 1991))
NOAA_YEARS = (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023, 2025)
BANDWIDTHS = np.array(
    [125.799765, 251.599530, 503.199060, 1006.398120, 2012.796241]
)
HISTORICAL_VARIABLES = {"d2m", "t2m", "sp"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def need(path: Path, purpose: str = "required artifact") -> Path:
    require(path.is_file(), f"Missing {purpose}: {path.relative_to(PROJECT)}")
    require(path.stat().st_size > 0, f"Empty {purpose}: {path.relative_to(PROJECT)}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, purpose: str = "audit") -> dict:
    return json.loads(need(path, purpose).read_text(encoding="utf-8"))


def read_csv(path: Path, purpose: str = "result", **kwargs) -> pd.DataFrame:
    return pd.read_csv(need(path, purpose), keep_default_na=False, **kwargs)


def verify_inventory(items: list[dict], label: str) -> None:
    require(items, f"{label} inventory is empty")
    for item in items:
        relative = item.get("file") or item.get("path")
        require(relative is not None, f"{label} inventory item lacks a path")
        path = need(PROJECT / relative, label)
        if "bytes" in item:
            require(path.stat().st_size == int(item["bytes"]),
                    f"Size mismatch in {label}: {relative}")
        if "sha256" in item:
            require(sha256(path) == item["sha256"],
                    f"Hash mismatch in {label}: {relative}")


def require_bandwidths(frame: pd.DataFrame, label: str) -> None:
    values = np.sort(frame["bandwidth_km"].astype(float).unique())
    require(len(values) == 5 and np.allclose(values, BANDWIDTHS, atol=1e-5, rtol=0),
            f"{label} does not contain the five frozen physical bandwidths")


def equirectangular_distance_km(lon1, lat1, lon2, lat2):
    mean_latitude = np.deg2rad((lat1 + lat2) / 2)
    dx = (lon1 - lon2) * 111.32 * np.cos(mean_latitude)
    dy = (lat1 - lat2) * 110.57
    return np.sqrt(dx * dx + dy * dy)


def verify_historical_acquisition_and_analysis() -> None:
    provenance = read_json(
        HIST_DATA / "historical_cds_point_provenance.json",
        "historical CDS provenance",
    )
    require(provenance.get("dataset") == "reanalysis-era5-land-timeseries",
            "Historical provenance names the wrong CDS product")
    require(tuple(provenance.get("years_retained", [])) == HISTORICAL_YEARS,
            "Historical provenance does not retain every summer from 1950 to 1990")
    require(provenance.get("raw_archives_retained") is True,
            "Historical raw archives were not marked as retained")
    archives = provenance.get("archives", [])
    require(len(archives) == 121, "Historical provenance must list 121 point archives")
    require(len({int(item["site_id"]) for item in archives}) == 121,
            "Historical archive site identifiers are not unique")
    for item in archives:
        path = need(PROJECT / item["file"], "historical CDS archive")
        require(zipfile.is_zipfile(path), f"Invalid historical ZIP: {path.name}")
        require(path.stat().st_size == int(item["bytes"]),
                f"Historical archive size mismatch: {path.name}")
        require(sha256(path) == item["sha256"],
                f"Historical archive hash mismatch: {path.name}")

    trimmed = sorted((HIST_DATA / "trimmed_points").glob("site_*_jja_buffers.nc"))
    require(len(trimmed) == 121, "Expected 121 historical trimmed point files")
    for path in trimmed:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 91512,
                    f"Unexpected historical point time dimension: {path.name}")
            require(set(dataset.data_vars) == HISTORICAL_VARIABLES,
                    f"Unexpected historical point variables: {path.name}")

    panels = sorted((HIST_DATA / "hourly_points").glob(
        "era5_land_*_jja_121sites.nc"
    ))
    require(len(panels) == 41, "Expected 41 historical yearly hourly panels")
    panel_years = []
    for path in panels:
        panel_years.append(int(path.name.split("_")[2]))
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 2232 and
                    dataset.sizes.get("site") == 121,
                    f"Unexpected historical panel dimensions: {path.name}")
            require(set(dataset.data_vars) == HISTORICAL_VARIABLES,
                    f"Unexpected historical panel variables: {path.name}")
    require(tuple(panel_years) == HISTORICAL_YEARS,
            "Historical hourly panels do not span exactly 1950--1990")

    daily_paths = sorted((HIST_DATA / "daily_fields").glob(
        "era5_land_*_jja_daily_fields.csv.gz"
    ))
    require(len(daily_paths) == 41, "Expected 41 historical daily-field files")
    daily_years = []
    for path in daily_paths:
        frame = read_csv(path, "historical daily field", usecols=[
            "year", "analysis_date", "day_definition", "site_id", "wbt"
        ])
        daily_years.append(int(frame.year.iloc[0]))
        require(len(frame) == 92 * 121,
                f"Historical daily-field row count is wrong: {path.name}")
        require(frame.year.nunique() == 1 and set(frame.day_definition) == {"utc"},
                f"Historical daily-field definition is wrong: {path.name}")
        require(frame.analysis_date.nunique() == 92 and frame.site_id.nunique() == 121,
                f"Historical daily field is incomplete: {path.name}")
        require(np.isfinite(frame.wbt.astype(float)).all(),
                f"Historical daily field has nonfinite WBT: {path.name}")
    require(tuple(daily_years) == HISTORICAL_YEARS,
            "Historical daily fields do not span exactly 1950--1990")

    build = read_json(HIST_OUT / "historical_field_build_manifest.json")
    analysis = read_json(HIST_OUT / "historical_analysis_manifest.json")
    require(tuple(build.get("years", [])) == HISTORICAL_YEARS,
            "Historical field-build manifest is incomplete")
    require(analysis.get("summers") == 41 and
            analysis.get("month_year_records") == 123 and
            analysis.get("daily_fields") == 3772 and
            analysis.get("sites") == 121,
            "Historical analysis manifest dimensions are wrong")
    require(analysis.get("randomisation", {}).get("draws") == 99999 and
            analysis.get("randomisation", {}).get("seed") == 20260810 and
            analysis.get("randomisation", {}).get(
                "specified_before_historical_value_retrieval") is True,
            "Historical product-shift design differs from the frozen protocol")
    climatology = analysis.get("climatology_reference", {})
    require(climatology.get("period") == "1950-1990" and
            climatology.get("modern_1991_2025_values_used") is False and
            "41 historical summers" in climatology.get("definition", ""),
            "Historical climatology reference is not period-specific")
    audit = read_csv(HIST_OUT / "historical_data_audit.csv")
    require(len(audit) == 41 and tuple(audit.year.astype(int)) == HISTORICAL_YEARS,
            "Historical data audit is incomplete")
    require((audit.sites.astype(int) == 121).all() and
            (audit.hours_with_buffer.astype(int) == 2232).all(),
            "Historical data audit dimensions are wrong")

    overall = read_csv(HIST_OUT / "historical_overall_results.csv")
    require(len(overall) == 1, "Historical overall result must have one row")
    row = overall.iloc[0]
    require(int(row.years) == 41 and int(row.negative_years) == 40,
            "Historical overall recurrence is not 40/41")
    require(abs(float(row.estimate) - (-0.110431)) < 1e-6 and
            float(row.ci_upper) < 0 and
            float(row.loo_min) < float(row.loo_max) < 0,
            "Historical overall estimate, interval, or leave-one-out range changed")

    scales = read_csv(HIST_OUT / "historical_scale_results.csv")
    require(len(scales) == 5, "Historical scale summary must have five rows")
    require_bandwidths(scales, "Historical scale summary")
    require((scales.years.astype(int) == 41).all() and
            list(scales.negative_years.astype(int)) == [38, 39, 40, 41, 41] and
            (scales.estimate.astype(float) < 0).all(),
            "Historical scale recurrence or direction changed")

    periods = read_csv(HIST_OUT / "historical_period_results.csv")
    require(list(periods.period) == ["1950-1978", "1979-1990"] and
            list(periods.years.astype(int)) == [29, 12] and
            list(periods.negative_years.astype(int)) == [28, 12] and
            (periods.ci_upper.astype(float) < 0).all(),
            "Historical forcing-period split is incomplete or changed")

    randomisation = read_csv(
        HIST_OUT / "historical_global_cyclic_randomisation.csv"
    )
    require(len(randomisation) == 1, "Historical global randomisation row is missing")
    rand = randomisation.iloc[0]
    require(int(rand.draws) == 99999 and int(rand.seed) == 20260810 and
            int(rand.random_values_at_or_below_observed) == 0 and
            abs(float(rand.p_value) - 1e-5) < 1e-12,
            "Historical global product-shift result changed")
    random_scales = read_csv(
        HIST_OUT / "historical_global_cyclic_randomisation_scales.csv"
    )
    require(len(random_scales) == 5, "Historical scale randomisation is incomplete")
    require_bandwidths(random_scales, "Historical scale randomisation")

    basis = read_csv(HIST_OUT / "historical_basis_decomposition.csv")
    require(len(basis) == 10 and
            set(basis.basis) == {"latitude", "latitude_longitude"},
            "Historical geographic-basis decomposition is incomplete")
    require_bandwidths(basis, "Historical basis decomposition")
    require((basis.structured_component_negative_years.astype(int) == 41).all() and
            basis.identity_error.astype(float).abs().max() <= 1e-12,
            "Historical basis recurrence or exact identity failed")
    energy = read_csv(HIST_OUT / "historical_energy_decomposition.csv")
    require(len(energy) == 5, "Historical energy decomposition is incomplete")
    require_bandwidths(energy, "Historical energy decomposition")
    require(energy.identity_error.astype(float).abs().max() <= 1e-12,
            "Historical energy decomposition identity failed")
    broad = energy.loc[np.isclose(energy.bandwidth_km.astype(float), BANDWIDTHS[-1])]
    require(len(broad) == 1 and
            int(broad.iloc[0].climatology_anomaly_cross_component_negative_years) == 41 and
            float(broad.iloc[0].anomaly_energy_component_ci_lower) < 0 <
            float(broad.iloc[0].anomaly_energy_component_ci_upper),
            "Historical broad-scale energy interpretation changed")
    print("[OK] historical acquisition, estimates, randomisation and decompositions")


def verify_continuous_calendar_and_ratio_stress() -> None:
    audit = read_json(METHOD_OUT / "extension_empirical_methods_audit.json")
    reproduction = audit.get("primary_reproduction", {})
    require(audit.get("site_count") == 121 and
            audit.get("summer_count_held_out") == 33,
            "Empirical-method audit dimensions are wrong")
    require(float(reproduction.get("absolute_error", 1)) <= 1e-12 and
            abs(float(reproduction.get("computed", 0)) - (-0.0728344831814084)) <= 1e-12,
            "Empirical-method extension does not reproduce the primary estimate")

    continuous = read_csv(METHOD_OUT / "extension_continuous_profile_summary.csv")
    require(len(continuous) == 1 and int(continuous.iloc[0].summers) == 33 and
            int(continuous.iloc[0].negative_summers) == 31 and
            float(continuous.iloc[0].estimate) < 0 and
            float(continuous.iloc[0].ci_upper) < 0,
            "Continuous log-Q profile result is incomplete or changed")
    continuous_scales = read_csv(
        METHOD_OUT / "extension_continuous_scale_summary.csv"
    )
    require(len(continuous_scales) == 5, "Continuous log-Q scale summary is incomplete")
    require_bandwidths(continuous_scales, "Continuous log-Q scale summary")
    require((continuous_scales.summers.astype(int) == 33).all() and
            list(continuous_scales.negative_summers.astype(int)) == [29, 30, 30, 32, 32] and
            (continuous_scales.ci_upper.astype(float) < 0).all(),
            "Continuous log-Q scale recurrence or intervals changed")
    months = read_csv(METHOD_OUT / "extension_continuous_month_slopes.csv")
    require(len(months) == 105 * 5 and months.record_id.nunique() == 105,
            "Continuous monthly slope panel is incomplete")
    require_bandwidths(months, "Continuous monthly slopes")

    coverage = read_csv(METHOD_OUT / "extension_calendar_matching_coverage.csv")
    require(list(coverage.window_days.astype(int)) == [3, 5] and
            (coverage.records.astype(int) == 99).all() and
            (coverage.high_days_total.astype(int) == 792).all() and
            list(coverage.high_days_matched.astype(int)) == [670, 750] and
            (coverage.records_with_no_match.astype(int) == 0).all(),
            "Calendar-match coverage counts changed")
    require(np.allclose(coverage.pooled_high_day_coverage.astype(float),
                        [670 / 792, 750 / 792], atol=1e-12, rtol=0),
            "Calendar-match coverage fractions do not match their counts")
    calendar = read_csv(METHOD_OUT / "extension_calendar_matching_profile_summary.csv")
    require(list(calendar.window_days.astype(int)) == [3, 5] and
            (calendar.summers.astype(int) == 33).all() and
            list(calendar.negative_summers.astype(int)) == [24, 25] and
            (calendar.ci_upper.astype(float) < 0).all(),
            "Calendar-match profile results are incomplete or changed")
    calendar_scales = read_csv(
        METHOD_OUT / "extension_calendar_matching_scale_summary.csv"
    )
    require(len(calendar_scales) == 10 and
            set(calendar_scales.window_days.astype(int)) == {3, 5},
            "Calendar-match scale summary is incomplete")
    require_bandwidths(calendar_scales, "Calendar-match scale summary")
    records = read_csv(METHOD_OUT / "extension_calendar_matching_records.csv")
    require(len(records) == 2 * 105 * 5 and
            records.record_has_no_match.astype(str).str.lower().eq("false").all(),
            "Calendar-match record panel is incomplete or contains empty matches")

    ratio_audit = read_json(METHOD_OUT / "extension_ratio_stress_audit.json")
    require(ratio_audit.get("seed") == 20260811 and
            ratio_audit.get("replications_per_cell") == 2000 and
            ratio_audit.get("product_shift_draws") == 999 and
            ratio_audit.get("paired_scenarios") is True,
            "Heavy-tail ratio stress design differs from the frozen protocol")
    ratio = read_csv(METHOD_OUT / "extension_ratio_stress_results.csv.gz")
    require(len(ratio) == 2 * 3 * 2000,
            "Heavy-tail ratio stress results do not contain 12,000 rows")
    require(set(ratio.scenario) == {"null", "alternative_minus_7pct"} and
            set(ratio.measure) == {"raw_ratio", "log_ratio", "bounded_symmetric"},
            "Heavy-tail stress scenarios or measures are incomplete")
    counts = ratio.groupby(["scenario", "measure"]).size()
    require((counts == 2000).all() and
            ratio.p_value.astype(float).between(0, 1).all(),
            "Heavy-tail stress cells or p-values are invalid")
    ratio_summary = read_csv(METHOD_OUT / "extension_ratio_stress_summary.csv")
    require(len(ratio_summary) == 6 and
            (ratio_summary.replications.astype(int) == 2000).all(),
            "Heavy-tail ratio summary is incomplete")
    for scenario in ("null", "alternative_minus_7pct"):
        cell = ratio_summary.loc[ratio_summary.scenario.eq(scenario)].set_index("measure")
        require(cell.loc["bounded_symmetric", "rmse"] <
                cell.loc["log_ratio", "rmse"] < cell.loc["raw_ratio", "rmse"],
                f"Heavy-tail RMSE ordering changed in {scenario}")
    null = ratio_summary.loc[ratio_summary.scenario.eq("null")]
    require((null.rejection_rate_05.astype(float) - 0.05).abs().max() <= 0.005,
            "Heavy-tail null randomisation size is no longer near 0.05")
    print("[OK] continuous intensity, date matching and heavy-tail stress test")


def verify_elevation_basis() -> None:
    audit = read_json(ELEV_OUT / "elevation_basis_audit.json")
    raw = need(PROJECT / audit["raw_file"], "ERA5-Land invariant geopotential")
    require(sha256(raw) == audit["raw_sha256"],
            "Invariant geopotential hash mismatch")
    require(audit.get("basis_spaces") == [
        "latitude", "latitude_longitude", "latitude_longitude_elevation"
    ], "Elevation audit basis spaces changed")
    require(float(audit.get("maximum_day_identity_error", 1)) <= 1e-12 and
            float(audit.get("maximum_summary_identity_error", 1)) <= 1e-12 and
            float(audit.get("maximum_primary_reproduction_error", 1)) <= 2e-10,
            "Elevation-basis identity or primary reproduction failed")
    verify_inventory(audit.get("outputs", []), "elevation-basis output")
    sites = read_csv(
        PROJECT / "data" / "era5_invariant" / "era5_land_121_site_elevation.csv"
    )
    require(len(sites) == 121 and sites.site_id.nunique() == 121,
            "Invariant elevation site table is incomplete")
    elevation_column = next((name for name in sites.columns if "elev" in name.lower()), None)
    require(elevation_column is not None and
            np.isfinite(sites[elevation_column].astype(float)).all(),
            "Invariant site elevation values are missing or nonfinite")

    summary = read_csv(ELEV_OUT / "elevation_basis_summary.csv")
    require(len(summary) == 15 and set(summary.basis) == {
        "latitude", "latitude_longitude", "latitude_longitude_elevation"
    }, "Elevation-basis summary does not have 3 x 5 cells")
    require_bandwidths(summary, "Elevation-basis summary")
    require((summary.total_effect_years.astype(int) == 33).all() and
            (summary.structured_component_years.astype(int) == 33).all() and
            summary.identity_error.astype(float).abs().max() <= 1e-12,
            "Elevation-basis summary years or identity failed")
    topo = summary.loc[summary.basis.eq("latitude_longitude_elevation")].sort_values(
        "bandwidth_km"
    )
    local, broad = topo.iloc[0], topo.iloc[-1]
    require(int(local.structured_component_negative_years) == 18 and
            float(local.structured_component_ci_lower) < 0 <
            float(local.structured_component_ci_upper),
            "Local topographic-basis result changed")
    require(int(broad.structured_component_negative_years) == 33 and
            float(broad.structured_component_ci_upper) < 0 and
            abs(float(broad.structured_component_estimate) - (-0.143361)) < 1e-6,
            "Broad topographic-basis result changed")
    yearly = read_csv(ELEV_OUT / "elevation_basis_year_components.csv")
    require(len(yearly) == 3 * 35 * 5 and yearly.year.nunique() == 35 and
            set(yearly.analysis_role) == {"development", "held_out"},
            "Elevation-basis yearly component panel is incomplete")
    held_out = yearly.loc[yearly.analysis_role.eq("held_out")]
    require(len(held_out) == 3 * 33 * 5 and held_out.year.nunique() == 33,
            "Elevation-basis held-out component panel is incomplete")
    print("[OK] invariant elevation and geographic-topographic basis sensitivity")


def verify_complete_simulation_tables() -> None:
    annual = read_csv(PROJECT / "output_corrected" /
                      "year_inference_simulation_summary.csv")
    keys = ["sample_size", "rho", "innovation", "effect"]
    require(len(annual) == 108 and not annual.duplicated(keys).any(),
            "Annual simulation summary does not contain 108 unique cells")
    require(set(annual.sample_size.astype(int)) == {20, 33, 60, 120} and
            set(annual.rho.astype(float)) == {0.0, 0.3, 0.6} and
            set(annual.innovation) == {"gaussian", "skewed", "t3"} and
            set(annual.effect) == {"null", "moderate", "application"},
            "Annual simulation factorial grid is incomplete")
    require((annual.simulations.astype(int) == 10000).all(),
            "Annual simulation cells do not all contain 10,000 replications")
    joint = read_csv(PROJECT / "output_corrected" /
                     "joint_dgp_simulation_summary.csv")
    expected_joint = {
        "joint_sign_null", "joint_gradient_contraction",
        "seasonal_progression_null", "seasonal_gradient_contraction",
        "anisotropic_null", "anisotropic_amplitude_contraction",
        "peak_selection_null", "peak_selection_contraction",
    }
    require(len(joint) == 8 and set(joint.scenario) == expected_joint and
            joint.scenario.nunique() == 8 and
            (joint.simulations.astype(int) == 1000).all(),
            "Coupled-DGP simulation summary is not the complete eight-cell design")

    audit = read_json(METHOD_OUT / "supp_complete_simulation_tables_audit.json")
    require(audit.get("year_cells") == 108 and audit.get("joint_cells") == 8,
            "Complete simulation-table audit has wrong cell counts")
    verify_inventory(audit.get("inputs", []), "complete-table input")
    verify_inventory(audit.get("outputs", []), "complete-table output")
    output_tex = need(METHOD_OUT / "supp_complete_simulation_tables.tex")
    portable_tex = need(PROJECT / "manuscript" / "generated" /
                        "supp_complete_simulation_tables.tex")
    require(sha256(output_tex) == sha256(portable_tex),
            "Portable and analysis copies of the complete simulation tables differ")
    text = portable_tex.read_text(encoding="utf-8")
    require(text.count(r"\begin{longtable}") == 3 and
            r"\label{tab:joint-dgp-complete}" in text and
            r"\label{tab:year-simulation-complete-a}" in text and
            r"\label{tab:year-simulation-complete-b}" in text and
            r"\label{tab:year-simulation-complete-c}" in text,
            "Generated supplement does not contain all four complete tables")
    print("[OK] complete 108-cell and eight-cell simulation tables")


def verify_cross_record_dependence_stress() -> None:
    summary = read_csv(
        METHOD_OUT / "extension_cross_record_stress_summary.csv",
        "cross-record dependence stress summary",
    )
    require(len(summary) == 4 and
            np.allclose(summary.shared_loading_lambda.astype(float),
                        [0.0, 0.3, 0.6, 0.8], atol=1e-12, rtol=0),
            "Cross-record stress test does not contain the four lambda cells")
    require((summary.replications.astype(int) == 2000).all() and
            (summary.product_shifts_per_replication.astype(int) == 999).all() and
            np.allclose(summary.nominal_alpha.astype(float), 0.05,
                        atol=1e-12, rtol=0),
            "Cross-record stress replication, shift, or alpha design changed")
    require(list(summary.rejections.astype(int)) == [103, 102, 99, 98] and
            np.allclose(summary.empirical_rejection_rate.astype(float),
                        [0.0515, 0.0510, 0.0495, 0.0490],
                        atol=1e-12, rtol=0),
            "Cross-record stress rejection results changed")
    expected_mcse = np.sqrt(
        summary.empirical_rejection_rate.astype(float) *
        (1 - summary.empirical_rejection_rate.astype(float)) / 2000
    )
    require(np.allclose(summary.monte_carlo_se.astype(float), expected_mcse,
                        atol=1e-12, rtol=0),
            "Cross-record stress Monte Carlo standard errors are inconsistent")
    correlations = summary.same_summer_july_august_daily_profile_correlation.astype(float)
    require(np.all(np.diff(correlations) > 0) and correlations.iloc[-1] > 0.75,
            "Cross-record shared-dependence diagnostic did not strengthen with lambda")

    audit = read_json(
        METHOD_OUT / "extension_cross_record_stress_audit.json",
        "cross-record dependence stress audit",
    )
    require(audit.get("seed") == 20260812 and
            audit.get("replications_per_lambda") == 2000 and
            audit.get("product_shift_draws_per_replication") == 999 and
            audit.get("nominal_alpha") == 0.05 and
            audit.get("site_count") == 121 and
            audit.get("shared_loadings_lambda") == [0.0, 0.3, 0.6, 0.8],
            "Cross-record stress audit design differs from the reported design")
    record_structure = audit.get("record_structure", {})
    require(record_structure.get("summers") == 2 and
            record_structure.get("records_per_summer") == 3 and
            record_structure.get("record_lengths") == [30, 31, 31, 30, 31, 31],
            "Cross-record stress record structure is wrong")
    checks = audit.get("numerical_checks", {})
    require(float(checks.get(
                "marginal_shared_process_one_step_rotation_max_abs_error", 1)) <= 1e-12 and
            float(checks.get(
                "unscaled_30_by_31_cross_covariance_relative_change_after_one_30_day_shift",
                0)) > 0.25 and
            float(checks.get(
                "unscaled_30_by_31_graph_dispersion_cross_covariance_relative_change_after_one_30_day_shift",
                0)) > 0.35,
            "Cross-record marginal invariance or joint-departure check failed")
    require("not an extension" in audit.get("analysis_role", "") and
            audit.get("paired_base_fields_labels_and_offsets_across_lambda") is True,
            "Cross-record stress role or paired design is not documented")
    verify_inventory(audit.get("inputs", []), "cross-record stress input")
    verify_inventory(audit.get("outputs", []), "cross-record stress output")
    need(PROJECT / "code" / "54_cross_record_dependence_stress.py",
         "cross-record dependence simulation script")
    print("[OK] cross-record dependence assumption stress test")


def verify_noaa_acquisition() -> tuple[pd.DataFrame, pd.DataFrame]:
    audit = read_json(
        NOAA_DATA / "noaa_extension_download_audit.json",
        "NOAA acquisition audit",
    )
    require(audit.get("history_snapshot_global_max_end") == 20250828 and
            audit.get("history_snapshot_in_scope_max_end") == 20250824 and
            audit.get("frozen_requested_end") == 20250831,
            "NOAA official-history snapshot boundaries changed")
    require(audit.get("original_candidate_count") == 0,
            "Frozen NOAA END>=2025-08-31 rule no longer records zero candidates")
    require(audit.get("operational_candidate_count") == 208 and
            audit.get("selected_stations") == 30 and
            audit.get("selection_uses_effects") is False,
            "NOAA administrative operationalisation or outcome-blind selection changed")
    require(tuple(audit.get("years", [])) == NOAA_YEARS and
            audit.get("minimum_complete_exact_hours") == 400 and
            audit.get("errors") == [],
            "NOAA acquisition years, threshold, or error log is invalid")
    operationalisation = audit.get("operationalization", "")
    require("2025-08-31" in operationalisation and
            "2025-08-24" in operationalisation and
            "zero candidates" in operationalisation and
            "outcome" in operationalisation.lower(),
            "NOAA acquisition audit does not explain the administrative-only change")
    history = need(PROJECT / audit["history_file"], "official ISD history snapshot")
    require(sha256(history) == audit["history_sha256"],
            "Official ISD history snapshot hash mismatch")
    history_frame = pd.read_csv(history)
    history_frame["BEGIN"] = pd.to_numeric(history_frame.BEGIN, errors="coerce")
    history_frame["END"] = pd.to_numeric(history_frame.END, errors="coerce")
    metadata_eligible = history_frame.loc[
        history_frame.LAT.between(20, 42) &
        history_frame.LON.between(105, 125) &
        history_frame.BEGIN.le(19920601) &
        history_frame.END.ge(20250824)
    ].copy()
    sites = pd.read_csv(PRIMARY_SITE_FILE)
    metadata_eligible["nearest_primary_site_km"] = [
        float(np.min(equirectangular_distance_km(
            row.LON, row.LAT, sites.lon, sites.lat
        )))
        for row in metadata_eligible.itertuples(index=False)
    ]
    retained_by_proximity = metadata_eligible.nearest_primary_site_km.le(150)
    require(len(metadata_eligible) == 228 and
            int(retained_by_proximity.sum()) == 208 and
            int((~retained_by_proximity).sum()) == 20 and
            abs(float(metadata_eligible.loc[
                retained_by_proximity, "nearest_primary_site_km"
            ].max()) - 147.73) < 0.01 and
            abs(float(metadata_eligible.loc[
                ~retained_by_proximity, "nearest_primary_site_km"
            ].min()) - 155.72) < 0.01,
            "NOAA 150-km proximity filter does not reproduce 228 -> 208")
    verify_inventory(audit.get("outputs", []), "NOAA acquisition output")

    manifest = read_csv(NOAA_DATA / "noaa_extension_station_manifest.csv",
                        dtype={"station": str, "USAF": str, "WBAN": str})
    require(len(manifest) == 30 and manifest.station.nunique() == 30 and
            list(manifest.selection_order.astype(int)) == list(range(1, 31)),
            "NOAA station manifest is not the frozen ordered 30-station set")
    require((manifest.administrative_end_rule.astype(int) == 20250824).all() and
            manifest.selection_rule.str.contains("no outcome", case=False).all(),
            "NOAA manifest does not preserve the outcome-blind administrative rule")
    require(manifest.LAT.astype(float).between(20, 42).all() and
            manifest.LON.astype(float).between(105, 125).all() and
            (manifest.nearest_analysis_site_km.astype(float) <= 150).all(),
            "NOAA manifest violates the frozen rectangle or 150-km support rule")

    qualification = read_csv(
        NOAA_DATA / "noaa_extension_station_year_qualification.csv",
        dtype={"station": str},
    )
    require(len(qualification) == 30 * 10 and
            not qualification.duplicated(["station", "year"]).any(),
            "NOAA station-year qualification table is incomplete or duplicated")
    require(set(qualification.station) == set(manifest.station) and
            set(qualification.year.astype(int)) == set(NOAA_YEARS),
            "NOAA qualification table does not cover the frozen stations and years")
    qualified_bool = qualification.qualified.astype(str).str.lower().eq("true")
    require((qualified_bool ==
             qualification.complete_exact_hours.astype(int).ge(400)).all(),
            "NOAA station-year qualification does not apply the 400-hour rule exactly")
    require((qualification.groupby("year")["station"].nunique() == 30).all() and
            qualified_bool.groupby(qualification.year).sum().min() >= 26,
            "NOAA station-year coverage is incomplete")
    for row in qualification.itertuples(index=False):
        source = need(PROJECT / row.source_file, "NOAA Global Hourly source file")
        require(sha256(source) == row.source_sha256,
                f"NOAA source hash mismatch: {row.source_file}")

    for year in NOAA_YEARS:
        panel = read_csv(
            NOAA_DATA / f"noaa_isd_extension_{year}_jja_exact_hours.csv.gz",
            "NOAA exact-hour evaluation panel",
            dtype={"station": str},
        )
        expected_rows = int(qualification.loc[
            qualification.year.astype(int) == year,
            "complete_exact_hours",
        ].sum())
        require(len(panel) == expected_rows and set(panel.year.astype(int)) == {year},
                f"NOAA {year} panel does not match its qualified complete-hour count")
        require(not panel.duplicated(["station", "time_utc"]).any(),
                f"NOAA {year} panel has duplicate station-hours")
        require(np.isfinite(panel[["temperature_c", "dewpoint_c",
                                   "station_pressure_pa"]].astype(float)).all().all(),
                f"NOAA {year} panel has nonfinite complete-hour values")
    print("[OK] NOAA acquisition and administrative operationalisation")
    return manifest, qualification


def verify_noaa_era_points(manifest: pd.DataFrame) -> set[str]:
    provenance = read_json(
        NOAA_ERA / "era5_station_point_provenance.json",
        "NOAA ERA5-Land station-point provenance (download is not optional)",
    )
    require(provenance.get("dataset") == "reanalysis-era5-land-timeseries" and
            tuple(provenance.get("years_retained", [])) == NOAA_YEARS,
            "NOAA ERA5-Land point provenance has wrong product or years")
    require(set(provenance.get("variables", [])) == {
        "2m_dewpoint_temperature", "2m_temperature", "surface_pressure"
    }, "NOAA ERA5-Land point provenance has wrong variables")
    require(provenance.get("station_manifest_sha256") == sha256(
        NOAA_DATA / "noaa_extension_station_manifest.csv"
    ), "NOAA ERA5-Land provenance uses a different station manifest")
    stations = provenance.get("stations", [])
    require(provenance.get("selected_stations") == 30 and
            provenance.get("available_station_series") == 28 and
            len(stations) == 30 and
            {str(item["station"]) for item in stations} == set(manifest.station),
            "NOAA ERA5-Land provenance does not contain the frozen 30 stations")
    available = {
        str(item["station"]) for item in stations if item["status"] == "available"
    }
    unavailable = [item for item in stations if item["status"] != "available"]
    require(len(available) == 28 and len(unavailable) == 2 and
            {str(item["station"]) for item in unavailable} == {
                "58472099999", "59562099999"
            } and
            all(item["status"] == "unavailable_nearest_land_mask_cell" and
                item["trimmed"] is None and "no finite" in item["reason"]
                for item in unavailable),
            "NOAA ERA5-Land land-mask exclusions changed")
    for item in stations:
        archive = need(PROJECT / item["archive"], "NOAA ERA5-Land raw archive")
        require(zipfile.is_zipfile(archive),
                f"Invalid NOAA ERA5-Land ZIP: {archive.name}")
        require(sha256(archive) == item["archive_sha256"],
                f"NOAA ERA5-Land archive hash mismatch: {archive.name}")
        if item["status"] != "available":
            continue
        trimmed = need(PROJECT / item["trimmed"], "NOAA ERA5-Land trimmed point")
        require(sha256(trimmed) == item["trimmed_sha256"],
                f"NOAA ERA5-Land trimmed hash mismatch: {trimmed.name}")
        with xr.open_dataset(trimmed, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 2232 * 10,
                    f"NOAA ERA5-Land trimmed time count is wrong: {trimmed.name}")
            require(set(dataset.data_vars) == HISTORICAL_VARIABLES,
                    f"NOAA ERA5-Land trimmed variables are wrong: {trimmed.name}")
            for variable in HISTORICAL_VARIABLES:
                require(bool(np.isfinite(dataset[variable].values).all()),
                        f"NOAA ERA5-Land trimmed values are nonfinite: {trimmed.name}")
    raw = sorted((NOAA_ERA / "raw_archives").glob("station_*.zip"))
    trimmed = sorted((NOAA_ERA / "trimmed_points").glob(
        "station_*_evaluation_jja.nc"
    ))
    require(len(raw) == 30 and len(trimmed) == 28,
            "NOAA ERA5-Land acquisition directories are incomplete")
    print("[OK] 28 finite NOAA-site ERA5-Land point series; two land-mask exclusions")
    return available


def verify_noaa_analysis(qualification: pd.DataFrame,
                         available: set[str]) -> None:
    audit = read_json(
        NOAA_OUT / "noaa_extension_analysis_audit.json",
        "NOAA extension analysis audit (analysis is not optional)",
    )
    require(tuple(audit.get("years", [])) == NOAA_YEARS and
            audit.get("qualified_station_year_threshold") == 400 and
            audit.get("minimum_stations_per_effect_field") == 10,
            "NOAA analysis protocol settings changed")
    require(np.allclose(sorted(audit.get("station_supported_bandwidths_km", [])),
                        BANDWIDTHS[2:], atol=1e-5, rtol=0),
            "NOAA analysis does not use the three frozen station-supported scales")
    require("frozen" in audit.get("effect_labels", "").lower(),
            "NOAA effect analysis does not identify the frozen ERA labels")
    require(audit.get("selected_stations") == 30 and
            audit.get("era_available_stations") == 28 and
            len(audit.get("era_unavailable_stations", [])) == 2 and
            audit.get("minimum_available_fields_per_regime_record") == 1,
            "NOAA analysis land-mask or record-availability rules changed")
    verify_inventory(audit.get("inputs", []), "NOAA analysis input")
    verify_inventory(audit.get("outputs", []), "NOAA analysis output")

    by_year = read_csv(NOAA_OUT / "noaa_extension_measurement_by_year.csv")
    overall = read_csv(NOAA_OUT / "noaa_extension_measurement_overall.csv")
    require(len(by_year) == 10 and tuple(by_year.year.astype(int)) == NOAA_YEARS and
            len(overall) == 1,
            "NOAA measurement summaries do not cover all ten years")
    metrics = ["bias_c", "mae_c", "rmse_c", "pooled_correlation",
               "within_station_centered_correlation",
               "equal_station_mean_correlation"]
    require(np.isfinite(by_year[metrics].astype(float)).all().all() and
            np.isfinite(overall[metrics].astype(float)).all().all(),
            "NOAA measurement agreement metrics are nonfinite")
    qualified_bool = qualification.qualified.astype(str).str.lower().eq("true")
    expected_hours = int(qualification.loc[
        qualified_bool & qualification.station.isin(available),
        "complete_exact_hours"
    ].astype(int).sum())
    require(int(overall.iloc[0].matched_hours) == expected_hours ==
            int(by_year.matched_hours.astype(int).sum()) == audit.get("matched_hours"),
            "NOAA matched-hour totals disagree with qualification or audit")

    matched = read_csv(NOAA_OUT / "noaa_extension_era5_matched.csv.gz",
                       dtype={"station": str})
    require(len(matched) == expected_hours and
            not matched.duplicated(["station", "time_utc"]).any(),
            "NOAA--ERA5-Land matched panel is incomplete or duplicated")
    require(np.isfinite(matched[["observed_wbt_c", "era_wbt_c"]].astype(float)).all().all(),
            "NOAA--ERA5-Land matched WBT values are nonfinite")

    fields = read_csv(NOAA_OUT / "noaa_extension_graph_fields.csv.gz")
    require_bandwidths(fields, "NOAA graph fields")
    require((fields.stations.astype(int) >= 10).all() and
            (fields.observed_q.astype(float) > 0).all() and
            (fields.era_q.astype(float) > 0).all() and
            set(fields.regime) <= {"high", "middle"},
            "NOAA graph fields violate support, positivity, or frozen regimes")
    events = read_csv(
        PROJECT / "output_confirmatory" / "sensitivity_event_manifest.csv"
    )
    events = events.loc[events.year.astype(int).isin(NOAA_YEARS)]
    field_times = set(pd.to_datetime(fields.time_utc))
    event_times = set(pd.to_datetime(events.peak_time))
    require(field_times <= event_times,
            "NOAA graph fields contain times outside the frozen ERA peak manifest")

    records = read_csv(NOAA_OUT / "noaa_extension_record_effects.csv")
    years = read_csv(NOAA_OUT / "noaa_extension_year_scale_effects.csv")
    scales = read_csv(NOAA_OUT / "noaa_extension_scale_summary.csv")
    profile = read_csv(NOAA_OUT / "noaa_extension_broad_profile_summary.csv")
    require(records.record_id.nunique() == 30 and len(records) == 30 * 5,
            "NOAA record-effect panel does not contain all 30 JJA records x five scales")
    require(len(years) == 10 * 5 and set(years.year.astype(int)) == set(NOAA_YEARS),
            "NOAA yearly-effect panel does not contain ten years x five scales")
    require(len(scales) == 5 and (scales.years.astype(int) == 10).all(),
            "NOAA scale summary does not contain five complete ten-year estimates")
    require_bandwidths(scales, "NOAA scale summary")
    require(len(profile) == 2 and set(profile.field) == {"NOAA station", "ERA5-Land"} and
            (profile.years.astype(int) == 10).all() and
            np.isfinite(profile[["estimate", "standard_error", "ci_lower",
                                 "ci_upper"]].astype(float)).all().all(),
            "NOAA broad-profile summary is incomplete or nonfinite")
    require(audit.get("effect_fields") == fields.time_utc.nunique() and
            audit.get("effect_records") == records.record_id.nunique(),
            "NOAA effect counts disagree with the analysis audit")
    print("[OK] NOAA--ERA5-Land measurement agreement and fixed-label effect analysis")


def verify_manuscript_integration() -> None:
    """Fail if completed extension analyses disappear from the submitted text."""
    main_text = need(MAIN_TEX, "main-manuscript source").read_text(encoding="utf-8")
    supplement = need(
        SUPPLEMENT_TEX, "supplement source"
    ).read_text(encoding="utf-8")

    required_main = [
        "The 1950--1990 extension gave",
        "The station comparison contains 175,172 exact-hour matches",
        "$-17.34\\%$",
        "both intervals include zero",
        "$\\pm3$-day estimator",
        "$\\pm5$-day estimator",
        "latitude--longitude--elevation basis",
        "2,000 heavy-tailed data sets",
        "figure07_noaa_agreement.pdf",
        "Supplementary Section~S11.4 studies one shared",
        "report its $p$-value as exploratory",
        "The protocol-defined three-component",
        "WGS84 geodesic distances",
        "leave-one-summer-out version gave",
        "fixed common support within each summer",
        "regional means over the available station network defined both",
        "station-based event definition on the available dynamic network",
    ]
    missing_main = [item for item in required_main if item not in main_text]
    require(not missing_main,
            f"Completed extension is missing from main manuscript: {missing_main}")
    stale_main = [
        "NOAA data assess two development summers",
        "The station check contains 35,512 matched hours",
    ]
    require(not any(item in main_text for item in stale_main),
            "Main manuscript retains the superseded two-summer NOAA account")

    required_supplement = [
        "\\section{Historical extension, 1950--1990}",
        "\\section{NOAA station extension in ten evaluation years}",
        "\\label{tab:noaa-extension-measurement}",
        "\\label{tab:noaa-extension-effect}",
        "175,172 paired hours",
        "\\section{Continuous intensity and calendar matching}",
        "\\section{Heavy-tailed ratio stress test}",
        "\\subsection{Cross-record dependence stress test}",
        "\\label{tab:cross-record-stress}",
        "primary finite-record summary",
        "Confirmatory component",
        "Confirmatory rule",
        "conditional product-invariance null",
        "Historical climatology--anomaly decomposition",
        "Thus $\\boldsymbol\\mu_m^{H}$ is the site-specific mean field",
        "leaving 208 candidates",
        "\\label{tab:domain-scale-curves}",
        "North edge inward $1.8^\\circ$, relabelled",
        "\\input{generated/supp_complete_simulation_tables.tex}",
    ]
    missing_supplement = [
        item for item in required_supplement if item not in supplement
    ]
    require(not missing_supplement,
            "Completed extension is missing from supplement: "
            f"{missing_supplement}")
    print("[OK] completed extension results are integrated in main and supplement")


def verify_post_review_sensitivities() -> None:
    """Verify the spatial-support, LOSO-climatology, and station checks."""
    output_dir = PROJECT / "output_revision_sensitivity"
    audit = read_json(
        output_dir / "revision_sensitivity_audit.json",
        "post-review sensitivity audit",
    )
    script = need(PROJECT / "code" / "55_revision_sensitivity_analyses.py")
    require(audit.get("script_sha256") == sha256(script),
            "Post-review sensitivity script hash changed")
    for item in audit.get("outputs", []):
        path = need(output_dir / item["file"], "post-review sensitivity output")
        require(sha256(path) == item["sha256"],
                f"Post-review sensitivity hash mismatch: {item['file']}")

    spatial = read_csv(output_dir / "revision_spatial_profile_summary.csv")
    require(len(spatial) == 26 and
            (spatial.n_years.astype(int) == 33).all() and
            (spatial.estimate.astype(float) < 0).all(),
            "Spatial-support sensitivity panel is incomplete or changes sign")
    curves = read_csv(output_dir / "revision_spatial_scale_curves.csv")
    require(len(curves) == 26 * 5 and
            curves.groupby("analysis").size().eq(5).all() and
            curves.groupby("analysis").estimate.apply(
                lambda values: values.astype(float).is_monotonic_decreasing
            ).all(),
            "Spatial-support curves do not retain all five ordered bandwidths")

    loso = read_csv(output_dir / "revision_climatology_loso_scale_curve.csv")
    require(len(loso) == 2 * 5 * 4 and
            float(audit["climatology_decomposition"]
                  ["leave_one_summer_out_monthly_climatology"]
                  ["maximum_daily_identity_error"]) <= 1e-12,
            "LOSO climatology decomposition is incomplete or fails identity")

    station = read_csv(output_dir / "revision_station_scale_uncertainty.csv")
    require(len(station) == 90 and
            station.groupby(["analysis", "field"]).size().eq(5).all(),
            "Station support/day-count/event-definition checks are incomplete")
    availability = read_csv(output_dir / "revision_station_availability.csv")
    all_years = availability.loc[availability.aggregation == "all_years"]
    require(len(all_years) == 2 and
            set(all_years.regime) == {"high", "middle"} and
            int(all_years.days_eligible_ge10.sum()) == 345,
            "State-specific station-availability audit changed")
    print("[OK] post-review spatial, climatology, and station sensitivities")


def main() -> None:
    verify_historical_acquisition_and_analysis()
    verify_continuous_calendar_and_ratio_stress()
    verify_elevation_basis()
    verify_complete_simulation_tables()
    verify_post_review_sensitivities()
    verify_cross_record_dependence_stress()
    manifest, qualification = verify_noaa_acquisition()
    available = verify_noaa_era_points(manifest)
    verify_noaa_analysis(qualification, available)
    verify_manuscript_integration()
    print("PASS: extension artifacts are present and internally consistent")


if __name__ == "__main__":
    main()
