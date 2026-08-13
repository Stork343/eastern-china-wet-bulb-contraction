#!/usr/bin/env python3
"""Verify the complete confirmatory data, analysis, and manuscript bundle."""

from __future__ import annotations

import hashlib
import json
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import xarray as xr
from pypdf import PdfReader


PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data" / "era5_confirmatory"
OUTPUT_DIR = PROJECT_DIR / "output_confirmatory"
PDF_FILE = PROJECT_DIR / "output" / "pdf" / \
    "humid_heat_spatial_contraction_jrssc.pdf"
SUPPLEMENT_PDF_FILE = PROJECT_DIR / "output" / "pdf" / \
    "humid_heat_spatial_contraction_jrssc_supplement.pdf"
MANUSCRIPT_FILE = PROJECT_DIR / "manuscript" / "main.tex"
SUPPLEMENT_FILE = PROJECT_DIR / "manuscript" / "supplement_theory.tex"
JRSSC_OUTPUT_DIR = PROJECT_DIR / "output_jrssc"
REPORT_FILE = OUTPUT_DIR / "completion_audit.json"
VARIABLES = {"d2m", "t2m", "sp", "u10", "v10", "swvl1", "ssrd"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def main() -> None:
    provenance_path = DATA_DIR / "cds_point_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    archives = provenance.get("archives", [])
    require(len(archives) == 121, "Provenance does not list 121 archives")
    for item in archives:
        path = PROJECT_DIR / item["file"]
        require(zipfile.is_zipfile(path), f"Invalid archive: {path.name}")
        require(path.stat().st_size == item["bytes"],
                f"Size mismatch: {path.name}")
        require(sha256(path) == item["sha256"],
                f"Hash mismatch: {path.name}")

    trimmed = sorted((DATA_DIR / "trimmed_points").glob(
        "site_*_jja_buffers.nc"
    ))
    require(len(trimmed) == 121, "Expected 121 trimmed point files")
    for path in trimmed:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 78120,
                    f"Unexpected time size: {path.name}")
            require(set(dataset.data_vars) == VARIABLES,
                    f"Unexpected variables: {path.name}")

    yearly = sorted((DATA_DIR / "hourly_points").glob(
        "era5_land_*_jja_121sites.nc"
    ))
    require(len(yearly) == 35, "Expected 35 yearly hourly panels")
    for path in yearly:
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            require(dataset.sizes.get("time") == 2232 and
                    dataset.sizes.get("site") == 121,
                    f"Unexpected panel dimensions: {path.name}")
            require(set(dataset.data_vars) == VARIABLES,
                    f"Unexpected panel variables: {path.name}")

    daily = sorted((DATA_DIR / "daily_fields").glob(
        "era5_land_*_jja_daily_fields.csv.gz"
    ))
    require(len(daily) == 35, "Expected 35 daily-field files")
    definitions = {"utc", "utc_plus_8", "utc_stull", "sitewise_max"}
    for path in daily:
        frame = pd.read_csv(path, usecols=[
            "year", "analysis_date", "day_definition", "site_id", "wbt"
        ])
        require(len(frame) == 4 * 92 * 121,
                f"Unexpected daily rows: {path.name}")
        require(set(frame.day_definition) == definitions,
                f"Unexpected definitions: {path.name}")
        counts = frame.groupby("day_definition").size()
        require((counts == 92 * 121).all(),
                f"Incomplete definition: {path.name}")
        require(frame.wbt.notna().all(), f"Missing WBT: {path.name}")

    audit = pd.read_csv(PROJECT_DIR / "output_corrected" /
                        "confirmatory_data_audit.csv")
    require(len(audit) == 35 and (audit.sites == 121).all() and
            (audit.hours_with_buffer == 2232).all(),
            "Confirmatory data audit is incomplete")

    primary = pd.read_csv(OUTPUT_DIR / "confirmatory_primary_results.csv")
    utc = primary.loc[primary.day_definition == "utc"]
    require(len(utc) == 1 and bool(utc.iloc[0].confirmatory_consistency),
            "UTC primary decision is missing or false")

    spatial_files = [
        "spatial_decomposition_role_maps.csv",
        "spatial_decomposition_identity_check.csv",
        "spatial_decomposition_scale_identity_check.csv",
    ]
    for name in spatial_files:
        require((OUTPUT_DIR / name).exists(),
                f"Spatial decomposition output is missing: {name}")
    identity = pd.read_csv(
        OUTPUT_DIR / "spatial_decomposition_identity_check.csv"
    )
    mapped_utc = identity.loc[identity.analysis_role == "confirmatory"]
    require(len(mapped_utc) == 1,
            "Confirmatory spatial identity check is missing")
    require(float(mapped_utc.iloc[0].absolute_error) <= 1e-10,
            "Mapped contributions do not recover the scalar effect")
    require(abs(float(mapped_utc.iloc[0].mapped_effect) -
                float(utc.iloc[0].estimate)) <= 1e-10,
            "Mapped and primary confirmatory effects disagree")
    spatial_maps = pd.read_csv(
        OUTPUT_DIR / "spatial_decomposition_role_maps.csv"
    )
    spatial_maps = spatial_maps.loc[
        spatial_maps.analysis_role == "confirmatory"
    ]
    scale_columns = [
        column for column in spatial_maps.columns
        if column.startswith("contribution_h_")
    ]
    scale_identity = pd.read_csv(
        OUTPUT_DIR / "spatial_decomposition_scale_identity_check.csv"
    )
    require(len(spatial_maps) == 121 and
            spatial_maps.site_id.nunique() == 121 and
            len(scale_columns) == 5,
            "Primary spatial maps do not retain 121 nodes and five scales")
    require(len(scale_identity) == 5 and
            scale_identity.absolute_error.max() <= 1e-10,
            "Scale-specific 121-node allocations fail their identities")
    require(max(abs(spatial_maps[scale_columns].sum().to_numpy() -
                    scale_identity.sort_values("scale_column").mapped_effect
                    .to_numpy())) <= 1e-10,
            "Retained scale-specific map columns disagree with their audit")

    smooth_surface_audit = pd.read_csv(
        JRSSC_OUTPUT_DIR / "primary_smooth_map_surface_audit.csv"
    )
    smooth_build_audit = json.loads((
        JRSSC_OUTPUT_DIR / "primary_smooth_map_build_audit.json"
    ).read_text())
    require(len(smooth_surface_audit) == 7 and
            (smooth_surface_audit.analysis_nodes == 121).all() and
            (smooth_surface_audit.prediction_cells > 10000).all() and
            (smooth_surface_audit.display_resolution_degrees == 0.1).all() and
            (smooth_surface_audit.spline_k == 30).all() and
            smooth_surface_audit.display_only.eq(True).all() and
            smooth_surface_audit.predictions_enter_estimation_or_inference
            .eq(False).all(),
            "Primary smooth maps are not audited as 121-node display-only fits")
    finite_smooth_identities = smooth_surface_audit.loc[
        smooth_surface_audit.reference_sum.notna(), "raw_identity_error"
    ]
    require(len(finite_smooth_identities) == 4 and
            finite_smooth_identities.max() <= 1e-10 and
            smooth_build_audit.get("prohibited_dense_inputs_used") is False and
            smooth_build_audit["display_interpolation"]
            ["predictions_enter_estimation_or_inference"] is False,
            "Primary map smoothing audit permits analytical interpolation")

    require(PDF_FILE.exists(), "Final manuscript PDF is missing")
    pdf = PdfReader(PDF_FILE)
    require(len(pdf.pages) <= 24,
            "Main manuscript exceeds the 24-page submission target")
    require(SUPPLEMENT_PDF_FILE.exists(),
            "Supplementary proof PDF is missing")
    supplement_pdf = PdfReader(SUPPLEMENT_PDF_FILE)

    manuscript = MANUSCRIPT_FILE.read_text()
    main_text = manuscript.split("\\begin{thebibliography}", maxsplit=1)[0]
    section_titles = re.findall(r"\\section\{([^}]+)\}", main_text)
    expected_sections = [
        "Introduction",
        "Multiscale graph regime contrasts",
        "Finite-record inference and study design",
        "Simulation study",
        "Humid-heat application",
        "Discussion",
    ]
    require(section_titles == expected_sections,
            "Main-text section order does not match the six-section design")
    abstract_match = re.search(
        r"\\abstract\{(.*?)\n\}\s*\n\s*\\keywords\{",
        manuscript,
        flags=re.DOTALL,
    )
    require(abstract_match is not None, "Abstract is missing")
    abstract_plain = re.sub(r"\\[A-Za-z]+(?:\[[^]]*\])?", " ",
                            abstract_match.group(1))
    abstract_words = re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*",
                                abstract_plain)
    require(len(abstract_words) <= 200,
            "Abstract exceeds the JRSS C 200-word limit")
    require("\\cite" not in abstract_match.group(1),
            "Abstract contains a citation command")
    require("product cyclic" not in abstract_match.group(1).lower() and
            "p_{" not in abstract_match.group(1) and
            "p=" not in abstract_match.group(1),
            "Abstract must not foreground the post-result product-shift p-value")
    require("weighted mean squared" in abstract_match.group(1) and
            "root-mean-square" in abstract_match.group(1) and
            "Gaussian graph bandwidth" in abstract_match.group(1) and
            "sparse external station comparison" in abstract_match.group(1),
            "Abstract lacks the squared/RMS distinction, bandwidth wording, or station limitation")
    abstract_abbreviations = [
        token for token in ("NOAA", "ERA5")
        if token in abstract_match.group(1)
    ]
    require(not abstract_abbreviations,
            f"Abstract contains abbreviations: {abstract_abbreviations}")
    require(
        "\\documentclass[namedate,webpdf,modern,mediumone]"
        "{oup-authoring-template}" in manuscript and
        "\\onecolumn" not in manuscript,
        "Official OUP modern medium single-column author-year template is missing",
    )
    require("\\graphicspath{{figures/}}" in manuscript and
            "../output_" not in manuscript,
            "Manuscript figure paths are not submission-portable")
    introduction = main_text.split("\\section{Introduction}", 1)[1].split(
        "\\section{Multiscale graph regime contrasts}", 1
    )[0]
    require("\\subsection" not in introduction,
            "Introduction must not contain subsections")
    require(main_text.count("\\begin{proposition}") == 0 and
            main_text.count("\\end{proposition}") == 0 and
            "the proof and the plus-one Monte Carlo" in main_text,
            "Technical product-shift proposition must remain in the supplement")
    require("\\begin{theorem}" not in main_text and
            "\\begin{proof}" not in main_text and
            "\\label{thm:asymptotic}" not in main_text and
            "\\appendix" not in manuscript,
            "The theorem, HAC details, and proofs must remain outside the main text")
    require(SUPPLEMENT_FILE.exists(),
            "Supplementary proof source is missing")
    supplement = SUPPLEMENT_FILE.read_text()
    require("primary finite-record summary" in supplement and
            "Pre-access extension" in supplement and
            "Exploratory conditional randomisation" in supplement and
            "conditional product-invariance null" in supplement and
            "Confirmatory component" in supplement and
            "Confirmatory rule" in supplement and
            "2.56\\times10^{-6}" in supplement and
            "3.54\\times10^{-6}" in supplement and
            "3.31\\times10^{-5}" in supplement and
            "Historical climatology--anomaly decomposition" in supplement and
            "process-level ratio argument additionally assumes" in supplement and
            "negative-moment condition in Equation" not in supplement and
            "Secondary reference summaries" not in supplement and
            "Extension plan set before retrieval" not in supplement,
            "Supplementary inferential roles or negative-moment wording are stale")
    require(supplement.count("\\begin{proposition}") == 1 and
            supplement.count("\\begin{theorem}") == 1 and
            supplement.count("\\begin{proof}") == 2 and
            supplement.count("\\end{proof}") == 2,
            "Supplement must contain the proposition, strong-mixing theorem, and two proofs")
    require("\\label{prop:product-shift}" in supplement and
            "\\label{thm:s-process}" in supplement and
            "\\label{eq:s-hac-estimator}" in supplement,
            "Supplementary proposition or process-level HAC reference is missing")
    proof_one = supplement.split(
        "\\section{Finite-record product cyclic-shift validity}", 1
    )[1].split("\\section{Process-level strong-mixing and HAC reference}", 1)[0]
    proof_two = supplement.split(
        "\\section{Process-level strong-mixing and HAC reference}", 1
    )[1].split("\\section{Global product-shift result}", 1)[0]
    require(len(re.findall(r"[A-Za-z]+", proof_one)) >= 650 and
            "eq:s-mc-exchangeable-rank" in proof_one and
            "p_h^{\\mathrm{adj}}" in proof_one,
            "Proposition proof lacks the orbit, Monte Carlo, or FWER argument")
    require(len(re.findall(r"[A-Za-z]+", proof_two)) >= 550 and
            "eq:s-davydov-bound" in proof_two and
            "eq:s-fixed-lag-limit" in proof_two,
            "Strong-mixing proof lacks covariance, HAC, or fixed-lag details")
    simulation_text = main_text.split("\\section{Simulation study}", 1)[1].split(
        "\\section{Humid-heat application}", 1
    )[0]
    application_text = main_text.split(
        "\\section{Humid-heat application}", 1
    )[1].split("\\section{Discussion}", 1)[0]
    require("report its $p$-value as exploratory" in main_text and
            "In the post-result exploratory calculation" in main_text and
            "The later product-shift result is exploratory" in main_text and
            "cyclic-shift result provides the finite-record inferential assessment"
            not in main_text,
            "Global product-shift timing or exploratory status is inconsistent")
    require(simulation_text.count("\\begin{table}") == 0 and
            application_text.count("\\begin{table}") == 1 and
            main_text.count("\\begin{table}") == 3 and
            "\\section*{Tables}" not in main_text,
            "Main tables must be contribution, status, and physical-scale tables")
    main_figure_count = main_text.count("\\begin{figure}")
    require(main_figure_count == 7,
            "Main text must contain exactly seven integrated figures")
    require(main_text.count("\\figalttext{") == main_figure_count,
            "Every main-text figure must use the OUP alt-text command")
    expected_main_figures = [
        "figure01_simulation_diagnostics.pdf",
        "figure02_study_area.pdf",
        "figure03_multiscale_evidence.pdf",
        "figure04_primary_spatial_decomposition.pdf",
        "figure05_energy_decomposition.pdf",
        "figure06_application_robustness.pdf",
        "figure07_noaa_agreement.pdf",
    ]
    included_main_figures = re.findall(
        r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", main_text
    )
    require(included_main_figures == expected_main_figures,
            "Main-text portable figure sequence is not the seven-figure design")
    for figure_block in re.findall(
            r"\\begin\{figure\}(.*?)\\end\{figure\}",
            main_text, flags=re.DOTALL):
        require(not re.search(r"\([a-d]\)", figure_block),
                "Figure panels must use uppercase labels for JRSS C")
    require("Maozai Tian\\ORCID{0000-0002-0515-4477}" in manuscript and
            "Corresponding author." in manuscript and
            "mztian@ruc.edu.cn" in manuscript,
            "Maozai Tian correspondence details or ORCID are missing")
    require("0000-0002-9248-6874" in manuscript and
            "0009-0001-1812-1834" in manuscript,
            "Jian Hou or Tan Meng ORCID is missing from the title page")

    sensitivity_required = [
        "sensitivity_robustness_summary.csv",
        "sensitivity_scale_diagnostics.csv",
        "sensitivity_temporal_diagnostics.csv",
        "sensitivity_gradient_decomposition.csv",
        "sensitivity_gradient_yearly.csv",
        "sensitivity_hourly_sensitivity.csv",
        "sensitivity_seasonal_timing.csv",
        "sensitivity_seasonal_timing_yearly.csv",
        "sensitivity_event_manifest.csv",
        "sensitivity_station_validation.csv",
        "sensitivity_station_spatial_validation.csv",
        "sensitivity_station_advanced_validation.csv",
        "sensitivity_station_graph_validation.csv",
    ]
    for name in sensitivity_required:
        require((OUTPUT_DIR / name).exists(),
                f"sensitivity diagnostic is missing: {name}")
    robustness = pd.read_csv(OUTPUT_DIR / "sensitivity_robustness_summary.csv")
    primary_review = robustness.loc[
        (robustness.analysis == "equal_site_ratio") &
        (robustness.estimand == "ratio_effect"),
        "estimate",
    ]
    require(len(primary_review) == 1 and
            abs(float(primary_review.iloc[0]) - float(utc.iloc[0].estimate))
            <= 1e-12,
            "sensitivity diagnostic does not reproduce the primary effect")
    scale_diag = pd.read_csv(OUTPUT_DIR / "sensitivity_scale_diagnostics.csv")
    require(len(scale_diag) == 5 and scale_diag["min"].min() > 2 and
            scale_diag.coefficient_of_variation.max() < 0.20,
            "Middle-day denominator diagnostics are unstable")
    gradient = pd.read_csv(
        OUTPUT_DIR / "sensitivity_gradient_decomposition.csv"
    )
    require(len(gradient) == 5 and gradient.identity_error.max() <= 1e-12,
            "Latitude-gradient decomposition fails its exact identity")
    require((gradient.gradient_negative_years == 33).all(),
            "Gradient recurrence is not reproduced across all held-out summers")
    hourly = pd.read_csv(OUTPUT_DIR / "sensitivity_hourly_sensitivity.csv")
    fixed_hour = hourly.loc[hourly.analysis == "fixed_hour"]
    require(len(fixed_hour) == 24 and (fixed_hour.estimate < 0).all() and
            (fixed_hour.ci_upper < 0).all(),
            "Fixed-hour sensitivity is incomplete or internally inconsistent")
    relabelled_hour = hourly.loc[
        hourly.analysis == "fixed_hour_relabelled"
    ]
    require(len(relabelled_hour) == 24 and
            (relabelled_hour.estimate < 0).all() and
            (relabelled_hour.ci_upper < 0).all(),
            "Fixed-hour relabelling is incomplete or internally inconsistent")
    seasonal = robustness.loc[
        robustness.analysis.isin([
            "site_year_month_linear_detrended",
            "leave_one_year_daily_climatology_anomaly",
        ]) & (robustness.estimand == "ratio_effect")
    ]
    require(len(seasonal) == 2 and (seasonal.estimate < 0).all() and
            (seasonal.ci_upper < 0).all(),
            "Seasonal-progression controls are incomplete")
    station_spatial = pd.read_csv(
        OUTPUT_DIR / "sensitivity_station_spatial_validation.csv"
    )
    require(station_spatial.pairs.sum() == 420216 and
            station_spatial.signed_pair_difference_correlation.min() > 0.95,
            "Station-pair spatial validation is incomplete")
    station_advanced = pd.read_csv(
        OUTPUT_DIR / "sensitivity_station_advanced_validation.csv"
    )
    residual_station = station_advanced.loc[
        station_advanced.diagnostic ==
        "distance_bin_residual_semivariance_correlation"
    ]
    require(len(residual_station) == 1 and
            float(residual_station.iloc[0].estimate) > 0.90 and
            int(residual_station.iloc[0].clusters) == 184,
            "Distance-adjusted station validation or date bootstrap is missing")
    station_graph = pd.read_csv(
        OUTPUT_DIR / "sensitivity_station_graph_validation.csv"
    )
    require((station_graph.analysis == "matched_hour_graph_dispersion").sum()
            == 5 and
            (station_graph.analysis == "development_event_contrast").sum()
            == 5,
            "Station graph-profile validation is incomplete")

    extended_required = [
        "extended_global_cyclic_randomisation.csv",
        "extended_global_cyclic_randomisation_scales.csv",
        "extended_scale_physical_summary.csv",
        "extended_effect_measure_summary.csv",
        "extended_energy_decomposition.csv",
        "extended_energy_decomposition_yearly.csv",
        "extended_basis_decomposition.csv",
        "extended_basis_decomposition_yearly.csv",
        "extended_dense_bandwidth_profile.csv",
    ]
    for name in extended_required:
        require((OUTPUT_DIR / name).exists(),
                f"extended analysis output is missing: {name}")

    sensitivity_dir = PROJECT_DIR / "output_revision_sensitivity"
    sensitivity_required = [
        "revision_sensitivity_summary.csv",
        "revision_spatial_scale_curves.csv",
        "revision_spatial_profile_summary.csv",
        "revision_domain_masks.csv",
        "revision_climatology_loso_scale_curve.csv",
        "revision_climatology_loso_profile_summary.csv",
        "revision_station_scale_uncertainty.csv",
        "revision_station_availability.csv",
        "revision_station_record_effects.csv",
        "revision_sensitivity_audit.json",
    ]
    for name in sensitivity_required:
        require((sensitivity_dir / name).exists(),
                f"Post-review sensitivity output is missing: {name}")

    sensitivity_audit = json.loads(
        (sensitivity_dir / "revision_sensitivity_audit.json").read_text()
    )
    require(
        sensitivity_audit.get("script_sha256") ==
        sha256(PROJECT_DIR / "code" / "55_revision_sensitivity_analyses.py"),
        "Post-review sensitivity script hash disagrees with its audit",
    )
    for item in sensitivity_audit.get("outputs", []):
        output_path = sensitivity_dir / item["file"]
        require(output_path.exists() and sha256(output_path) == item["sha256"],
                f"Post-review sensitivity hash mismatch: {item['file']}")

    spatial_sensitivity = pd.read_csv(
        sensitivity_dir / "revision_spatial_profile_summary.csv"
    ).set_index("analysis")
    spatial_expected = {
        "primary_equirect_equal_fixed_labels": -0.0728344831814079,
        "primary_wgs84_equal_fixed_labels": -0.06973712319708796,
        "primary_wgs84_area_fixed_labels": -0.07802418189600338,
        "primary_wgs84_area_relabelled": -0.07369593119050659,
        "natural_earth_china_mainland_intersection_wgs84_equal_domain_mean_relabelled":
            -0.0635290360129743,
    }
    for analysis, expected in spatial_expected.items():
        require(analysis in spatial_sensitivity.index and
                abs(float(spatial_sensitivity.loc[analysis, "estimate"]) -
                    expected) <= 1e-12,
                f"Post-review spatial result changed: {analysis}")
    boundary_relabelled = spatial_sensitivity.loc[
        spatial_sensitivity.label_rule == "domain_mean_relabelled"
    ]
    require(boundary_relabelled.estimate.max() < 0 and
            boundary_relabelled.estimate.min() >= -0.08 and
            boundary_relabelled.estimate.max() <= -0.063,
            "Boundary relabelling no longer preserves the reported range")

    loso = pd.read_csv(
        sensitivity_dir / "revision_climatology_loso_scale_curve.csv"
    )
    loso_broad = loso.loc[
        (loso.analysis == "leave_one_summer_out_monthly_climatology") &
        (loso.bandwidth_km > 2000)
    ].set_index("component")
    require(abs(float(loso_broad.loc[
                "anomaly_energy_component", "estimate"]) -
                0.003988234249113046) <= 1e-12 and
            abs(float(loso_broad.loc[
                "climatology_anomaly_cross_component", "estimate"]) +
                0.13667231408658762) <= 1e-12,
            "Leave-one-summer-out climatology decomposition changed")

    station_sensitivity = pd.read_csv(
        sensitivity_dir / "revision_sensitivity_summary.csv"
    )
    station_sensitivity = station_sensitivity.loc[
        (station_sensitivity.section == "station_comparison") &
        (station_sensitivity.field_or_component == "observed_wbt_c")
    ].set_index("analysis")
    station_expected = {
        "dynamic_support_frozen_era_labels_min1": -0.17368174971664208,
        "dynamic_support_frozen_era_labels_min2": -0.1342232326809351,
        "dynamic_support_frozen_era_labels_min3": -0.15244151163448613,
        "dynamic_support_frozen_era_labels_min5": -0.1455542563510763,
        "dynamic_support_station_defined_peak_and_labels": -0.15988788778923385,
        "year_specific_fixed_common_support_frozen_era_labels": -0.18076925378583827,
    }
    for analysis, expected in station_expected.items():
        require(analysis in station_sensitivity.index and
                abs(float(station_sensitivity.loc[analysis, "estimate"]) -
                    expected) <= 1e-12,
                f"Post-review station result changed: {analysis}")

    global_shift = pd.read_csv(
        OUTPUT_DIR / "extended_global_cyclic_randomisation.csv"
    )
    require(len(global_shift) == 1,
            "Global product-shift output must contain one joint test")
    global_row = global_shift.iloc[0]
    require(int(global_row.draws) == 99999 and
            int(global_row.month_year_records) == 99 and
            int(global_row.bandwidths) == 5 and
            int(global_row.random_values_at_or_below_observed) == 0,
            "Global product-shift design is not the 99,999-draw joint test")
    require(abs(float(global_row.p_value) - 1e-5) <= 1e-15 and
            abs(float(global_row.minimum_attainable_p) - 1e-5) <= 1e-15 and
            abs(float(global_row.observed) - float(utc.iloc[0].estimate))
            <= 1e-12,
            "Global p-value or observed five-scale statistic is inconsistent")
    global_scales = pd.read_csv(
        OUTPUT_DIR / "extended_global_cyclic_randomisation_scales.csv"
    )
    require(len(global_scales) == 5 and
            global_scales.notna().all().all() and
            global_scales.bandwidth_km.is_monotonic_increasing,
            "Scale-wise global randomisation output is incomplete")

    physical = pd.read_csv(
        OUTPUT_DIR / "extended_scale_physical_summary.csv"
    )
    require(len(physical) == 5 and physical.notna().all().all(),
            "Physical Q/RMS output must contain five complete rows with no NA")
    rms_middle_error = (
        physical.rms_middle_c.pow(2) - 2 * physical.q_middle_c2
    ).abs().max()
    rms_high_error = (
        physical.rms_high_c.pow(2) - 2 * physical.q_high_c2
    ).abs().max()
    require(max(rms_middle_error, rms_high_error) <= 1e-10,
            "Reported physical RMS values do not equal sqrt(2Q)")

    energy = pd.read_csv(
        OUTPUT_DIR / "extended_energy_decomposition.csv"
    )
    energy_sum = (
        energy.anomaly_energy_component_estimate +
        energy.climatology_anomaly_cross_component_estimate +
        energy.climatology_energy_component_estimate
    )
    require(len(energy) == 5 and energy.notna().all().all() and
            (energy.total_effect_estimate - energy_sum).abs().max() <= 1e-12 and
            energy.identity_error.max() <= 1e-12 and
            energy.climatology_energy_component_estimate.abs().max() <= 1e-15,
            "Five-scale climatology-anomaly energy identity fails")
    energy_yearly = pd.read_csv(
        OUTPUT_DIR / "extended_energy_decomposition_yearly.csv"
    )
    energy_yearly_sum = (
        energy_yearly.anomaly_energy_component +
        energy_yearly.climatology_anomaly_cross_component +
        energy_yearly.climatology_energy_component
    )
    require(len(energy_yearly) == 35 * 5 and
            energy_yearly.notna().all().all() and
            (energy_yearly.analysis_role == "held_out").sum() == 33 * 5 and
            (energy_yearly.analysis_role == "development").sum() == 2 * 5 and
            (energy_yearly.total_effect - energy_yearly_sum).abs().max()
            <= 1e-12,
            "Annual climatology-anomaly energy identity fails")

    basis = pd.read_csv(OUTPUT_DIR / "extended_basis_decomposition.csv")
    basis_yearly = pd.read_csv(
        OUTPUT_DIR / "extended_basis_decomposition_yearly.csv"
    )
    require(len(basis) == 2 * 5 and basis.notna().all().all() and
            basis.identity_error.max() <= 1e-12 and
            len(basis_yearly) == 2 * 35 * 5 and
            basis_yearly.notna().all().all() and
            (basis_yearly.analysis_role == "held_out").sum() == 2 * 33 * 5 and
            (basis_yearly.analysis_role == "development").sum() == 2 * 2 * 5 and
            (basis_yearly.total_effect - basis_yearly.structured_component -
             basis_yearly.residual_component).abs().max() <= 1e-12,
            "Latitude/planar basis allocations fail their exact identities")

    bandwidth_profile = pd.read_csv(
        OUTPUT_DIR / "extended_dense_bandwidth_profile.csv"
    )
    bandwidth_curve = bandwidth_profile.loc[
        bandwidth_profile.source == "31_point_curve"
    ].sort_values("bandwidth_km")
    require(len(bandwidth_curve) == 31 and
            bandwidth_curve.bandwidth_km.nunique() == 31 and
            bandwidth_curve.notna().all().all() and
            bandwidth_curve.bandwidth_km.is_monotonic_increasing and
            bandwidth_curve.estimate.is_monotonic_decreasing,
            "Dense bandwidth result is not the complete monotone 31-point curve")

    convergence = pd.read_csv(
        PROJECT_DIR / "output_dense" / "extended_spatial_convergence.csv"
    )
    convergence_summary = pd.read_csv(
        PROJECT_DIR / "output_dense" /
        "extended_spatial_convergence_summary.csv"
    )
    expected_configurations = {
        "primary_121": 121,
        "longitude_refined": 239,
        "latitude_refined": 236,
        "dense_465": 465,
    }
    observed_configurations = dict(
        convergence[["configuration", "sites"]].drop_duplicates().itertuples(
            index=False, name=None
        )
    )
    require(len(convergence) == 4 * 5 and
            convergence.notna().all().all() and
            observed_configurations == expected_configurations and
            (convergence.groupby("configuration").size() == 5).all() and
            len(convergence_summary) == 4 and
            convergence_summary.notna().all().all() and
            set(convergence_summary.configuration) ==
            set(expected_configurations) and
            convergence_summary.loc[
                convergence_summary.configuration == "dense_465",
                ["max_absolute_difference", "rmse"],
            ].abs().max().max() <= 1e-15,
            "Four-configuration spatial convergence output is incomplete")

    # Preserve the literal design label ``null``; pandas otherwise treats it
    # as a missing-value token and silently breaks the null-cell audit.
    graph_sim = pd.read_csv(
        PROJECT_DIR / "output_corrected" / "graph_simulation_summary.csv",
        keep_default_na=False,
    )
    require(len(graph_sim) == 25 and
            (graph_sim.simulations == 1000).all() and
            (graph_sim.family == "alternative").sum() == 15 and
            "rejection_variogram_profile" in graph_sim.columns,
            "Spatial simulation output is not the 25-cell, 1,000-run design")
    joint_sim = pd.read_csv(
        PROJECT_DIR / "output_corrected" / "joint_dgp_simulation_summary.csv",
        keep_default_na=False,
    )
    require(len(joint_sim) == 8 and
            (joint_sim.simulations == 1000).all() and
            (joint_sim.family == "null").sum() == 4,
            "Joint-DGP simulation output is not the 8-cell, 1,000-run design")
    irregular_sim = pd.read_csv(
        PROJECT_DIR / "output_corrected" /
        "irregular_year_simulation_summary.csv"
    )
    require(len(irregular_sim) == 6 and
            (irregular_sim.simulations == 10000).all() and
            irregular_sim.bias.abs().max() < 0.001,
            "Irregular-calendar simulation is incomplete or biased")
    year_sim = pd.read_csv(
        PROJECT_DIR / "output_corrected" /
        "year_inference_simulation_summary.csv"
    )
    require(len(year_sim) == 108 and
            (year_sim.simulations == 10000).all(),
            "Year-level simulation output is not the 108-cell, 10,000-run design")
    dense_scale_identity = pd.read_csv(
        PROJECT_DIR / "output_dense" / "dense_scale_identity_check.csv"
    )
    dense_identity = pd.read_csv(
        PROJECT_DIR / "output_dense" / "dense_spatial_identity_check.csv"
    )
    require(len(dense_scale_identity) == 5 and
            dense_scale_identity.absolute_error.max() <= 1e-10 and
            dense_identity.absolute_error.max() <= 1e-10,
            "Dense-grid spatial allocations fail their exact identities")

    portable_figures = {
        JRSSC_OUTPUT_DIR / "fig4_simulation_diagnostics.pdf":
            "figure01_simulation_diagnostics.pdf",
        PROJECT_DIR / "output_corrected" / "fig5_study_area.pdf":
            "figure02_study_area.pdf",
        JRSSC_OUTPUT_DIR / "fig3_multiscale_evidence.pdf":
            "figure03_multiscale_evidence.pdf",
        JRSSC_OUTPUT_DIR / "fig_primary_spatial_decomposition.pdf":
            "figure04_primary_spatial_decomposition.pdf",
        JRSSC_OUTPUT_DIR / "fig_energy_decomposition.pdf":
            "figure05_energy_decomposition.pdf",
        JRSSC_OUTPUT_DIR / "fig5_application_robustness.pdf":
            "figure06_application_robustness.pdf",
        JRSSC_OUTPUT_DIR / "fig6_noaa_agreement.pdf":
            "figure07_noaa_agreement.pdf",
        JRSSC_OUTPUT_DIR / "supp_dense_bandwidth_profile.pdf":
            "supp_dense_bandwidth_profile.pdf",
        JRSSC_OUTPUT_DIR / "supp_spatial_convergence.pdf":
            "supp_spatial_convergence.pdf",
    }
    for generated_path, submission_name in portable_figures.items():
        require(generated_path.exists(),
                f"Generated figure is missing: {generated_path.name}")
        portable_path = (PROJECT_DIR / "manuscript" / "figures" /
                         submission_name)
        require(
            portable_path.exists(),
            f"Portable manuscript figure is missing: {submission_name}",
        )
        if submission_name in {
                "figure03_multiscale_evidence.pdf",
                "figure04_primary_spatial_decomposition.pdf"}:
            require(sha256(generated_path) == sha256(portable_path),
                    f"Portable smooth map differs from {generated_path.name}")
    supplement_figures = re.findall(
        r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", supplement
    )
    require(supplement_figures == [
        "supp_dense_bandwidth_profile.pdf",
        "supp_spatial_convergence.pdf",
    ], "Supplementary figure sequence is not the two-figure design")

    report = {
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "archives": len(archives),
        "archive_bytes": sum(item["bytes"] for item in archives),
        "trimmed_point_files": len(trimmed),
        "yearly_hourly_panels": len(yearly),
        "daily_field_files": len(daily),
        "data_years": [1991, 2025],
        "confirmatory_years": int(utc.iloc[0].years),
        "primary_profile_effect": float(utc.iloc[0].estimate),
        "primary_consistency": True,
        "spatial_identity_error": float(mapped_utc.iloc[0].absolute_error),
        "spatial_scale_identity_error": float(
            scale_identity.absolute_error.max()
        ),
        "spatial_map_sites": 121,
        "spatial_map_display_cells": int(
            smooth_surface_audit.prediction_cells.iloc[0]
        ),
        "spatial_map_spline_k": int(smooth_surface_audit.spline_k.iloc[0]),
        "spatial_map_predictions_enter_inference": False,
        "manuscript_pages": len(pdf.pages),
        "supplement_pages": len(supplement_pdf.pages),
        "manuscript_sections": section_titles,
        "manuscript_abstract_words": len(abstract_words),
        "main_text_figures": main_figure_count,
        "main_text_tables": main_text.count("\\begin{table}"),
        "global_product_shift_draws": int(global_row.draws),
        "global_product_shift_p_value": float(global_row.p_value),
        "energy_identity_error": float(energy.identity_error.max()),
        "dense_bandwidth_points": len(bandwidth_curve),
        "spatial_convergence_configurations": len(expected_configurations),
        "physical_rms_identity_error": float(max(
            rms_middle_error, rms_high_error
        )),
        "sensitivity_diagnostics": len(sensitivity_required),
        "sensitivity_denominator_minimum": float(scale_diag["min"].min()),
        "sensitivity_gradient_identity_error": float(
            gradient.identity_error.max()
        ),
        "station_spatial_pairs": int(station_spatial.pairs.sum()),
        "spatial_simulation_datasets": int(graph_sim.simulations.sum()),
        "joint_dgp_simulation_datasets": int(joint_sim.simulations.sum()),
        "irregular_year_simulation_datasets": int(
            irregular_sim.simulations.sum()
        ),
        "year_simulation_replications": int(year_sim.simulations.sum()),
        "dense_scale_identity_error": float(
            dense_scale_identity.absolute_error.max()
        ),
        "post_review_sensitivity_outputs": len(sensitivity_required),
        "wgs84_primary_profile": float(spatial_sensitivity.loc[
            "primary_wgs84_equal_fixed_labels", "estimate"
        ]),
        "china_land_profile": float(spatial_sensitivity.loc[
            "natural_earth_china_mainland_intersection_wgs84_equal_domain_mean_relabelled",
            "estimate",
        ]),
        "loso_broad_cross_component": float(loso_broad.loc[
            "climatology_anomaly_cross_component", "estimate"
        ]),
        "station_fixed_support_profile": float(station_sensitivity.loc[
            "year_specific_fixed_common_support_frozen_era_labels", "estimate"
        ]),
        "manuscript_sha256": sha256(PDF_FILE),
        "supplement_sha256": sha256(SUPPLEMENT_PDF_FILE),
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
