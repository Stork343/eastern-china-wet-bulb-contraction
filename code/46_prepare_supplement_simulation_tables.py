#!/usr/bin/env python3
"""Create portable supplementary tables for every retained simulation cell."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parent.parent
YEAR_FILE = PROJECT / "output_corrected" / "year_inference_simulation_summary.csv"
JOINT_FILE = PROJECT / "output_corrected" / "joint_dgp_simulation_summary.csv"
OUTPUT = PROJECT / "output_extension_methods"
GENERATED = PROJECT / "manuscript" / "generated"
OUTPUT_TEX = OUTPUT / "supp_complete_simulation_tables.tex"
PORTABLE_TEX = GENERATED / "supp_complete_simulation_tables.tex"
AUDIT_FILE = OUTPUT / "supp_complete_simulation_tables_audit.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pct(value: float, digits: int = 1) -> str:
    return f"{100 * float(value):.{digits}f}"


def make_year_tables(frame: pd.DataFrame) -> list[str]:
    expected = 4 * 3 * 3 * 3
    if len(frame) != expected:
        raise ValueError(f"Expected {expected} annual cells, found {len(frame)}")
    keys = ["sample_size", "rho", "innovation", "effect"]
    if frame.duplicated(keys).any():
        raise ValueError("Annual simulation cells are not unique")
    if set(frame.effect) != {"null", "moderate", "application"}:
        raise ValueError(f"Unexpected effect labels: {sorted(set(frame.effect))}")
    frame["effect"] = pd.Categorical(
        frame.effect, categories=["null", "moderate", "application"], ordered=True
    )
    frame = frame.sort_values(keys).reset_index(drop=True)

    lines = [
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\begin{longtable}{rrllrrrrr}",
        r"\caption{Complete 108-cell repeated-summer results: estimation and coverage. Bias and RMSE are percentage points; coverage is per cent.}\label{tab:year-simulation-complete-a}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & Bias & RMSE & Student & Lag 2 & Growing \\",
        r" & & & & \multicolumn{2}{c}{percentage points} & \multicolumn{3}{c}{coverage (\%)} \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{9}{l}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & Bias & RMSE & Student & Lag 2 & Growing \\",
        r"\midrule",
        r"\endhead",
        r"\midrule\multicolumn{9}{r}{Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for row in frame.itertuples(index=False):
        innovation = {"gaussian": "Gaussian", "skewed": "Skewed", "t3": "$t_3$"}[row.innovation]
        effect = {"null": "Null", "moderate": r"$-3.5\%$", "application": r"$-7\%$"}[row.effect]
        lines.append(
            f"{row.sample_size:d} & {row.rho:.1f} & {innovation} & {effect} & "
            f"{100 * row.bias:.3f} & {100 * row.rmse:.3f} & "
            f"{pct(row.coverage_student)} & {pct(row.coverage_nw2)} & "
            f"{pct(row.coverage_hac)} \\\\"
        )
    lines += [r"\end{longtable}", ""]

    lines += [
        r"\scriptsize",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{longtable}{rrllrrrr}",
        r"\caption{Complete 108-cell repeated-summer results: standardised-mean diagnostics. $L$ is the growing Bartlett lag; the standardised columns use the known long-run variance.}\label{tab:year-simulation-complete-b}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & $L$ & $z$ mean & $z$ SD & $z_{.025},z_{.975}$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{8}{l}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & $L$ & $z$ mean & $z$ SD & $z_{.025},z_{.975}$ \\",
        r"\midrule",
        r"\endhead",
        r"\midrule\multicolumn{8}{r}{Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for row in frame.itertuples(index=False):
        innovation = {"gaussian": "Gaussian", "skewed": "Skewed", "t3": "$t_3$"}[row.innovation]
        effect = {"null": "Null", "moderate": r"$-3.5\%$", "application": r"$-7\%$"}[row.effect]
        lines.append(
            f"{row.sample_size:d} & {row.rho:.1f} & {innovation} & {effect} & "
            f"{row.growing_lag:d} & {row.z_mean:.3f} & {row.z_sd:.3f} & "
            f"{row.z_q025:.2f}, {row.z_q975:.2f} \\\\"
        )
    lines += [r"\end{longtable}", ""]

    lines += [
        r"\scriptsize",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{longtable}{rrllrrrrr}",
        r"\caption{Complete 108-cell repeated-summer results: one-sided rejection percentages. The three-test column requires Student, lag-2 and sign values all to be at most 0.025. Under non-null cells the entries are power, not size.}\label{tab:year-simulation-complete-c}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & Student & Lag 2 & Growing & Sign & Three-test \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{9}{l}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"$R$ & $\rho$ & Innovation & Effect & Student & Lag 2 & Growing & Sign & Three-test \\",
        r"\midrule",
        r"\endhead",
        r"\midrule\multicolumn{9}{r}{Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for row in frame.itertuples(index=False):
        innovation = {"gaussian": "Gaussian", "skewed": "Skewed", "t3": "$t_3$"}[row.innovation]
        effect = {"null": "Null", "moderate": r"$-3.5\%$", "application": r"$-7\%$"}[row.effect]
        lines.append(
            f"{row.sample_size:d} & {row.rho:.1f} & {innovation} & {effect} & "
            f"{pct(row.rejection_student)} & {pct(row.rejection_nw2)} & "
            f"{pct(row.rejection_hac)} & {pct(row.rejection_sign)} & "
            f"{pct(row.rejection_three_025)} \\\\"
        )
    lines += [r"\end{longtable}", ""]
    return lines


def make_joint_table(frame: pd.DataFrame) -> list[str]:
    if len(frame) != 8 or frame.scenario.nunique() != 8:
        raise ValueError("Joint-DGP summary must contain eight unique scenarios")
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\caption{Complete results for the eight coupled label--field stress designs. Targets, bias and RMSE are percentage points; rejection is per cent over 1,000 replications. Graph and five-bin variogram profiles use identical fields, labels and product shifts.}",
        r"\label{tab:joint-dgp-complete}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3.2pt}",
        r"\begin{tabular}{P{0.19\textwidth}P{0.13\textwidth}rrrrrrrr}",
        r"\toprule",
        r"Scenario & Family & \multicolumn{4}{c}{Graph profile} & \multicolumn{4}{c}{Five-bin variogram} \\",
        r" & & Target & Bias & RMSE & Reject & Target & Bias & RMSE & Reject \\",
        r"\midrule",
    ]
    labels = {
        "joint_sign_null": "Shared-latent null",
        "joint_gradient_contraction": "Shared-latent contraction",
        "seasonal_progression_null": "Seasonal-progression null",
        "seasonal_gradient_contraction": "Seasonal-gradient contraction",
        "anisotropic_null": "Anisotropic null",
        "anisotropic_amplitude_contraction": "Anisotropic amplitude contraction",
        "peak_selection_null": "Peak-selection null",
        "peak_selection_contraction": "Peak-selection contraction",
    }
    for row in frame.itertuples(index=False):
        family = "Null" if row.family == "null" else "Alternative"
        lines.append(
            f"{labels[row.scenario]} & {family} & {100 * row.target_graph:.2f} & "
            f"{100 * row.bias_graph:.2f} & {100 * row.rmse_graph:.2f} & "
            f"{pct(row.rejection_graph)} & {100 * row.target_variogram:.2f} & "
            f"{100 * row.bias_variogram:.2f} & {100 * row.rmse_variogram:.2f} & "
            f"{pct(row.rejection_variogram)} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return lines


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    GENERATED.mkdir(parents=True, exist_ok=True)
    year = pd.read_csv(YEAR_FILE, keep_default_na=False)
    joint = pd.read_csv(JOINT_FILE, keep_default_na=False)
    text = "\n".join([
        "% Generated by code/46_prepare_supplement_simulation_tables.py",
        *make_joint_table(joint),
        *make_year_tables(year),
    ]) + "\n"
    OUTPUT_TEX.write_text(text, encoding="utf-8")
    PORTABLE_TEX.write_text(text, encoding="utf-8")
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(PROJECT)),
        "year_cells": len(year),
        "joint_cells": len(joint),
        "inputs": [
            {"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
            for path in (YEAR_FILE, JOINT_FILE)
        ],
        "outputs": [
            {"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
            for path in (OUTPUT_TEX, PORTABLE_TEX)
        ],
    }
    AUDIT_FILE.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(f"Generated {PORTABLE_TEX} with {len(year)} annual and {len(joint)} joint cells")


if __name__ == "__main__":
    main()
