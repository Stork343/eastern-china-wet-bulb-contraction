# Reproducibility materials — September application version

This archive contains the manuscript sources and PDFs, executable analysis
programs, protocols, spatial manifests and retained numerical results used
in the September 2026 application revision. The August public release at
https://github.com/Stork343/eastern-china-wet-bulb-contraction/releases/tag/jrssc-submission-v1-2026-08-13
is an earlier version. The current companion archive is
`jrssc_reproducibility_bundle.zip`, available with the matching manuscripts
and portable LaTeX source at
https://github.com/Stork343/eastern-china-wet-bulb-contraction/releases/tag/jrssc-submission-v5-2026-09-08.

## Inspecting the reported results

The root `output_revision_application/` contains:

- `config.json` and `revision_protocol.md`: the fixed application specification.
- `cell_weights.csv`, `adjacency.csv` and `day_manifest.parquet`: spatial
  support, connectivity, synchronous times and original event labels.
- `daily_footprints.parquet`, `regime_footprints.csv`,
  `annual_footprints.csv` and `coverage_summary.csv`: complete threshold curves
  and the high/middle comparisons, including the 24°C numerical example.
- `fold_manifest.csv`, `heldout_predictions.parquet`,
  `annual_comparison.csv` and `comparison_summary.csv`: the shared temporal
  splits, predictions and errors behind the model-comparison table.
- `paired_improvement_summary.csv`, `threshold_paired_improvement.csv`,
  `coverage_sensitivity.csv` and `threshold_grid_sensitivity.csv`: paired
  comparisons and the specified sensitivity checks.
- `software_versions.json`: package versions used for the application analysis.

The main table averages all 77 thresholds from 9.5 to 28.5°C, equally within
months and then over 33 evaluation summers. All-day errors cover 3,036 daily
fields; high-day errors cover 792 original upper-quartile days. Development
years 2015 and 2022 are excluded. The represented land area is
3,220,939.0300219203 km². Differences are calculated before rounding.

`manuscript/generated/` contains the numerical tables and application text
included in the manuscript. `submission/main_manuscript.pdf` and
`submission/supplementary_material.pdf` are the matching PDFs.

## Recomputing the application comparison

Run instructions are in `code/revision_202609/README.md`. That directory
contains `common.py`, the identical versioned protocol, and all seven entry
points from `00_audit.py` through `06_figures_tables.py`. The README explains
how to preserve the supplied outputs before a fresh recomputation.

Required constructed inputs are the original files under
`data/era5_confirmatory/daily_fields/` and `data/era5_dense/daily_fields/`.
Provider downloads and these large reconstructed panels are excluded from the
archive. Acquisition and construction programs are supplied as:

- `code/download_confirmatory_cds_points.py`
- `code/build_confirmatory_fields.py`
- `code/download_dense_cds_points.py`
- `code/build_dense_daily_fields.py`

The dense builder retains both `primary_grid_peak` and `dense_grid_peak`;
the September analysis explicitly selects `primary_grid_peak`.
Copernicus access requires the user's own credentials and acceptance of the
provider terms. Provider requests and file hashes are recorded in the included
provenance files and `output_revision_application/input_hashes.json`.
The Natural Earth land boundary is supplied by R package spData; its hash is
recorded in the fixed configuration. The Python dependencies are NumPy,
pandas, SciPy, scikit-learn, pyarrow, geopandas, shapely, pyproj and matplotlib.

## Earlier analyses retained with this version

The original numerical outputs are under `results/`; the primary, dense-grid
and historical-extension protocols remain at the archive root. Original
primary and dense-grid protocol copies are retained in `protocols/archive/`;
the root copies update only script filenames. The primary protocol hash cited
in the Supplement identifies the original archived copy. Their analysis
programs are under `code/`, including `confirmatory_analysis.R`,
`dense_resolution_analysis.R`, `analyze_historical_extension.py`,
`analyze_noaa_extension.py`, `extended_analyses.py` and
`revision_sensitivity_analyses.py`. Existing simulations, historical results
and station comparisons are retained without rerunning them for the
September application comparison. Legacy verifier scripts concern those
original analyses and should not be treated as checks of the current
manuscript page counts or added sections.

The portable LaTeX archive `jrssc_source_package.zip` includes all figures and
all generated TeX inputs. It can be compiled independently of the data tree.
