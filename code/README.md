# Corrected analysis pipeline

The manuscript uses only the scripts listed below. Run them from the repository
root in this order:

| Order | Script | Output |
|---:|---|---|
| 1 | `02_rebuild_consistent_era5.R` | Six pressure-aware WBT files on one fixed 121-site grid in `data/era5_consistent/` |
| 2 | `06_corrected_empirical.R` | Primary estimates, month-year effects, shift tests, and continuous model in `output_corrected/` |
| 3 | `04_corrected_null_calibration.R` | Null calibration of the circular-shift test |
| 4 | `07_corrected_sensitivity.R` | Threshold, daily-field, WBT-algorithm, and kernel sensitivity analyses |
| 5 | `13_corrected_figures.R` | Supplementary null-calibration figure |
| 6 | `16_graph_esh_empirical.R` | Five-bandwidth graph-dispersion profile, omnibus test, and FWER-adjusted scale tests |
| 7 | `17_verify_graph_theory.R` | Numerical verification of the graph identities and spectral bounds |
| 8 | `18_graph_esh_simulation.R` | Cyclic and AR(1) calibration plus three spatial mechanisms, comparing the graph profile with a five-bin variogram profile and five classical summaries |
| 9 | `20_define_spatial_grid.R` | Reconstruct and audit the 169-candidate/121-land-site sampling rule |
| 10 | `21_download_confirmatory_cds_points.py` | Acquire and assemble the frozen 1991-2025 ERA5-Land point panel from the official CDS time-series service |
| 11 | `22_build_confirmatory_fields.py` | Apply quality gates, pressure-aware WBT, and UTC/UTC+8 daily-field definitions |
| 12 | `23_confirmatory_analysis.R` | Run the held-out multi-year analysis with year as the replication unit |
| 13 | `24_year_level_inference_simulation.R` | Stress-test and summarise the 33-year confirmatory decision rule |
| 14 | `25_study_area_figure.R` | Generate the audited study-area and sampling-lattice figure |
| 15 | `26_download_noaa_isd.py` | Download a deterministic, dispersed NOAA ISD station sample |
| 16 | `27_validate_noaa_isd.R` | Match quality-controlled station hours to native ERA5-Land fields |
| 17 | `31_spatial_field_decomposition.R` | Decompose the profile and each graph scale exactly over the 121 sites, with raw-node identity checks |
| 18 | `39_sensitivity_diagnostics.py` | Recompute denominator, transformation, fixed and relabelled area weights, matched kernels, seasonal progression, fixed-hour relabelling, gradient uncertainty, and distance-adjusted station diagnostics |
| 19 | `40_joint_dgp_simulation.R` | Stress-test the graph and five-bin variogram profiles under shared latent weather, within-month progression, anisotropic covariance and peak-hour selection |
| 20 | `41_extended_analyses.py` | Run the 33-summer global product shift, physical-effect summaries, exact climatology--anomaly and nested-basis decompositions, dense bandwidth curve, and four-level spatial convergence |
| 21 | `38_jrssc_main_figures.R` | Assemble the unified vector figures and invoke the 121-node display renderer after all simulations and diagnostics finish |
| 22 | `30_verify_confirmatory_bundle.py` | Verify archives, panels, diagnostics, simulation scales, portable source paths, manuscript structure, figures, spatial identities, and final PDFs |

After the manuscript review, run `55_revision_sensitivity_analyses.py` before
the final verifiers. It writes the WGS84 and boundary-support curves, the
leave-one-summer-out climatology decomposition, and the station availability,
fixed-support, stricter-day-count, station-event and pressure checks to
`output_revision_sensitivity/`. These are explicitly additional sensitivity
analyses and do not replace the protocol-defined primary target.
At the end of the figure build, script 38 sources
`53_primary_smooth_spatial_surfaces.R`. Script 53 fits low-rank thin-plate REML
surfaces to the exact 121-node anomaly and allocation summaries, overlays all
source nodes, and writes the final Figures 3 and 4. The 0.1-degree raster is a
display grid only. It does not enter graph construction, estimation,
uncertainty calculations or node-sum identities. Script 53 can also be run
directly after script 31 when only these two figures need to be rebuilt.

## Nested-grid spatial-resolution sensitivity

The 465-site branch is secondary and does not replace the frozen 121-site
confirmatory analysis. Its grid contains every primary site exactly and adds
344 ERA5-Land land locations.

| Order | Script | Output |
|---:|---|---|
| 1 | `32_define_dense_spatial_grid.R` | Freeze and audit the nested 0.8 by 0.9 degree grid |
| 2 | `33_download_dense_cds_points.py` | Reuse 121 primary archives and acquire only the 344 new core-variable point series |
| 3 | `34_build_dense_daily_fields.py` | Build fields at primary-grid and dense-grid UTC peak times |
| 4 | `35_dense_resolution_analysis.R` | Reproduce the 121-site result, then run fixed-label and recomputed-label 465-site analyses at the original absolute bandwidths |
| 5 | `36_dense_spatial_surfaces.R` | Verify the dense node decomposition and generate clipped thin-plate REML display surfaces with sampled sites overlaid |
| 6 | `37_verify_dense_bundle.py` | Verify the nested manifest, new archives, dense panels, daily fields, embedded-primary reproduction, spatial identities, and figures |

After both branches and the sensitivity diagnostics finish, run
`41_extended_analyses.py`, then
`38_jrssc_main_figures.R` to assemble the unified main-text PDF figures in
`output_jrssc/`, then compile the manuscript and run
`30_verify_confirmatory_bundle.py`.

The smooth surfaces are descriptive. All estimates and tests use observed
node values; the profile and scale-specific contribution sums are checked
against the scalar graph estimates before plotting. Large record-level map and
simulation-draw tables are not retained because the fixed scripts reproduce
them from the archived data.

`esh_utils.R` contains the shared data checks, spatial metrics, classification,
resampling, and record-level summaries.
`graph_esh_utils.R` contains graph construction, normalized Laplacian
dispersion, profile effects, and product cyclic-shift inference.

## Definitions used by the manuscript

- One observation is the synchronous 121-site field at the hour of each UTC
  day's maximum regional-mean WBT.
- Event labels use mean WBT only and are computed within each month-year:
  upper quartile = high, interquartile range = middle.
- Spatial variance is the complete-graph baseline. The main extension uses
  normalized graph dispersion at 126, 252, 503, 1,006, and 2,013 km.
- The graph-profile summary is the equal-month, equal-scale and equal-summer
  mean relative high-minus-middle effect over those five prespecified
  bandwidths. Scale-specific adjusted values use
  the minimum relative effect from each simultaneous product shift.
- Development-record calculations use 9,999 independent within-record
  circular shifts. The primary 33-summer finite-record test uses 99,999
  product shifts across all 99 held-out month-year records while keeping each
  five-scale daily profile intact. Exact randomisation validity requires cyclic
  invariance of the joint outcome sequence conditional on the mean-WBT record.
  The linear AR(1) null is a calibration boundary case, not a proof of validity
  under stationarity.
- Primary WBT is pressure-aware and calculated by solving the saturated Bolton
  equivalent-potential-temperature relation. The Stull (2011) equation is kept
  as a sensitivity check.

## Confirmatory extension

The frozen decisions for the multi-year extension are in
`../CONFIRMATORY_ANALYSIS_PROTOCOL.md`. Years 2015 and 2022 are development
data; the other 33 years in 1991--2025 form the held-out sample. The main
paper reports the prespecified finite-record mean and its scale profile. The
Student, fixed-lag and sign calculations in the historical protocol are kept
as an audit trail. The manuscript now reports all three values and that their
prespecified consistency rule passed, but does not present the rule as
definitive long-run process inference. The protocol's stated freeze date and
SHA-256 are retained; because the first visible public commit is later, public
Git history alone is not treated as proof of the freeze time. The subsequently
added global product shift is labelled exploratory throughout.

Run `20_define_spatial_grid.R` before remote acquisition. The CDS point
downloader requires the packages in `requirements-era5.txt` and the personal
access token in `~/.cdsapirc`. It retains each raw point archive, resumes from
completed sites, trims the prespecified seasonal buffers, and assembles the 35
yearly panels consumed by script 22. The service normally permits only one
queued time-series request per account, so the default `--workers 1` should be
kept unless the CDS limit changes.

## Post-analysis 1950--1990 temporal extension

The frozen design is in `../EXTENSION_ANALYSIS_PROTOCOL.md`. This branch is
isolated from `data/era5_confirmatory/` and `output_confirmatory/`; it neither
rewrites nor reclassifies the 1991--2025 analysis.

| Order | Script | Output |
|---:|---|---|
| 1 | `42_download_historical_extension_cds_points.py` | Core-variable archives, trimmed JJA buffers and 41 yearly 121-site panels in `data/era5_historical_extension/` |
| 2 | `43_build_historical_extension_fields.py` | Pressure-aware Bolton WBT at the frozen UTC regional-peak hour and file-hash audit |
| 3 | `44_analyze_historical_extension.py` | Overall, five-scale, annual and forcing-segment effects; the 99,999-draw product shift; energy and latitude/planar basis decompositions in `output_historical_extension/` |

The downloader uses the same fixed site manifest but only the three variables
needed for WBT. Completed ZIP and NetCDF files are resumable and written
atomically. For a full run:

```sh
python3 code/42_download_historical_extension_cds_points.py --workers 4
python3 code/43_build_historical_extension_fields.py
python3 code/44_analyze_historical_extension.py
```

ECMWF documents that 1950--1978 ERA5-Land was forced by the preliminary ERA5
back extension, including a sub-optimal representation of some tropical
cyclones. Script 44 therefore reports 1950--1978 and 1979--1990 separately;
the split is fixed by production history rather than the observed effect.

The completed branch contains 121 retained point archives, 41 yearly panels,
123 month-year records and 3,772 daily fields. Its equal-summer estimate is
-11.0431% (40/41 summers negative), and the five-scale profile runs from
-5.9566% to -17.1735%. The 99,999-draw product shift uses seed 20260810 and
returns the plus-one value `p=0.00001`. The retained numerical results and
hash manifests are under `output_historical_extension/`. The dated pre-access
protocol explicitly includes that product shift. The historical energy
decomposition estimates its monthly climatology only from the 41 summers in
1950--1990 and is reported as an exploratory structural analysis because it
was not listed in the pre-access extension protocol.

## Post-analysis estimator and structural checks

These scripts use the frozen 1991--2025 fields and do not change the primary
event labels or estimand.

| Order | Script | Output |
|---:|---|---|
| 1 | `44_extension_empirical_methods.py` | Continuous log-energy slopes and same-month +/-3-day and +/-5-day matching in `output_extension_methods/` |
| 2 | `45_ratio_stress_test.R` | Paired raw-ratio, log-ratio and bounded-contrast stress tests under the null and -7% alternative |
| 3 | `48_elevation_basis_extension.py` | Official invariant-geopotential acquisition and latitude, latitude-longitude and latitude-longitude-elevation basis decompositions in `output_elevation_basis/` |
| 4 | `46_prepare_supplement_simulation_tables.py` | Complete 108-cell repeated-summer and eight-cell joint-DGP LaTeX tables plus a source/output hash audit |
| 5 | `54_cross_record_dependence_stress.py` | Targeted size stress test for shared phase-aligned temporal structure that preserves marginal record-level cyclic stationarity but violates joint product invariance |

The continuous profile slope is -0.065047 per degree C (31/33 summers
negative). The +/-3-day and +/-5-day effects are -5.2108% and -5.7237%, with
84.596% and 94.697% high-day coverage; neither window leaves a whole record
unmatched. The ratio experiment uses 2,000 replications per scenario, 999
product shifts and seed 20260811. Under the null, raw, log and bounded
rejection rates are 4.70%, 4.90% and 4.75%, while their RMSE values are
1.4725, 0.3283 and 0.2771.

The cross-record experiment uses the actual 30/31/31 JJA record lengths for
two summers, 2,000 paired replications at shared loadings 0, 0.3, 0.6 and 0.8,
and 999 product shifts per test (seed 20260812). The respective lower-tail
rejection rates are 5.15%, 5.10%, 4.95% and 4.90%. Its covariance audit checks
record-level rotation invariance numerically and confirms that an independent
month shift changes both the shared 30-by-31-day field cross-covariance and the
corresponding graph-dispersion cross-covariance when the loading is positive.
This is a sensitivity calculation, not an extension of the product-invariance
proposition.

The invariant geopotential file is retained at
`data/era5_invariant/era5_land_geopotential_20200101.nc`. At 2,013 km the
latitude-longitude-elevation structured component is -14.3361% and is negative
in 33/33 summers; at 126 km it is -1.1641% and its reference interval crosses
zero. The calculation treats the three-column space jointly and does not
assign an order-specific elevation effect.

The complete simulation table is written both to
`output_extension_methods/supp_complete_simulation_tables.tex` and the
portable manuscript path
`manuscript/generated/supp_complete_simulation_tables.tex`. Its audit requires
exactly 108 repeated-summer cells and eight joint-DGP cells.

## Non-development-year NOAA extension

The station branch uses the enumerated set
`{1992+4k: k=0,...,7} union {2023,2025}` and excludes the 2015 and 2022
development years. The archived protocol records the set but no further
outcome-derived sampling rule.

| Order | Script | Output |
|---:|---|---|
| 1 | `47_download_noaa_extension.py` | Official history snapshot, outcome-blind 30-station maximin manifest, ten yearly exact-hour panels and station-year qualification audit under `data/noaa_isd_extension/` |
| 2 | `49_download_noaa_extension_era5_points.py` | Matching ERA5-Land point archives, 28 finite trimmed series and provenance for the two unavailable land-mask cells |
| 3 | `50_analyze_noaa_extension.py` | Measurement agreement and frozen-label graph effects in `output_noaa_extension/` |

The official history snapshot did not reach the frozen `END >= 20250831`
threshold: its in-scope maximum was 20250824, so the literal rule produced no
candidates. Script 47 records both counts and uses the in-scope snapshot
maximum only for this administrative end field. The evaluation years,
rectangle, 150-km rule, station count and outcome-blind maximin selection are
unchanged. The 150-km rule is the minimum equirectangular distance to a primary
121-site node and is applied before maximin selection. It reduces 228
metadata-eligible stations to 208 candidates and then 30 selected stations;
`selection_uses_effects` is `false` in the audit. The cutoff amendment and
station selection precede retrieval of station WBT observations, matched
ERA5-Land series and graph contrasts.

ERA5-Land has no finite nearest-cell values for Shengsi (`58472099999`) and
Taidong (`59562099999`). Their archives and reasons are retained; neither was
replaced after effect inspection. The final 28-station panel contains 175,172
matched hours. Its bias, MAE, RMSE and within-station-centred correlation are
-0.2626 degrees C, 0.9383 degrees C, 1.2703 degrees C and 0.9160. The equal
503/1,006/2,013-km effect is -17.3355% for NOAA (9/10 years negative) and
-20.9859% for matched ERA5-Land (10/10 negative).

## Extension verification

Run the extension verifier after scripts 42--50:

```sh
python3 code/51_verify_extension_bundle.py
```

It checks the historical archives and estimates, continuous and calendar
methods, ratio-stress and cross-record replication counts, invariant elevation,
all 108+8 simulation cells, NOAA administrative operationalisation, both land-mask
exclusions, measurement summaries and frozen-label effects. The design
definitions and analysis-plan hash are recorded in
`../EXTENSION_ANALYSIS_PROTOCOL.md` and the
retained audit JSON files.

## Submission archives

After the two PDFs and submission-facing documents are final, run:

```sh
python3 code/52_build_submission_packages.py
```

Script 52 writes the portable LaTeX source archive and the scientific
reproducibility archive under `submission/`. The source archive contains only the two
TeX sources, the generated 108+8 simulation table, nine vector figures and its
README. The scientific archive adds code, protocols, retained outputs and small
provider manifests while excluding the large raw ERA5-Land archives and
reconstructed panels. Both archives are written atomically and checked with
`ZipFile.testzip()` before replacement.
