# Reproducibility materials

This repository contains the material needed to inspect the retained results and
reconstruct the analysis: executable scripts, dated protocols, grid and station
manifests, small provider-provenance files, processed results, simulation
summaries, manuscript sources, final figures and the two submitted PDFs. It
excludes submission correspondence, author administration, private working
notes, obsolete duplicates, machine-specific files and the approximately
11-GB local provider-data tree.

The large ERA5-Land and NOAA source files remain with their official
providers. The download and construction scripts rebuild them in the paths
expected by the analysis scripts. The retained provenance files record the
requests, station selection and source-file hashes.

## Repository access

The scientific archive is larger than the JRSS C limit of 2 MB for an
individual supplementary-material file and should not be uploaded as a single
supplementary file. Its submission-matched contents are frozen in the release
<https://github.com/Stork343/eastern-china-wet-bulb-contraction/releases/tag/jrssc-submission-v1-2026-08-13>.

## Directory guide

- `code/` contains the acquisition, construction, analysis, figure and
  verification programs. Numeric prefixes give the execution order.
- `data/` contains spatial grids, station manifests, compact station extracts
  and provenance needed to reconstruct provider data.
- `output_*` directories contain the retained numerical results and audit
  summaries used by the manuscript and supplement.
- `manuscript/` contains portable LaTeX sources, the generated complete
  simulation table and all submitted vector figures.
- `output/pdf/` contains the final main and supplementary PDFs.
- `SUBMISSION_MANIFEST.sha256` hashes the protocols, submitted text and PDFs,
  key result tables, verification programs and added sensitivity audit.
- The three protocol files at the archive root record the primary, dense-grid
  and extension designs.

Public filenames describe the analytical task. In particular,
`39_sensitivity_diagnostics.py` contains the transformation, weighting,
seasonal, hourly and station diagnostics, while `41_extended_analyses.py`
contains the global product shift, physical-scale summaries, exact energy and
basis decompositions, dense bandwidth profile and spatial-support comparison.
The additional checks requested before submission are kept separately in
`55_revision_sensitivity_analyses.py` and `output_revision_sensitivity/` so
their timing cannot be confused with the protocol-defined primary analysis.
They cover WGS84 distances and boundary support, leave-one-summer-out
climatology, station availability and fixed support, station-defined events,
and pressure conversion.

The Version 2 confirmatory protocol reports a freeze date of 2 August 2026 and
has SHA-256
`13c3a32532921e9230d00e88e07bb1aaa95e1a27ac970d843739bbb3090a3a70`.
The first commit visible in the public Git history is later, so that history
does not independently establish the stated freeze date. The manuscript and
supplement report the original three-component rule and all three results,
while clearly identifying the global product-shift analysis as exploratory.

## Analysis order

The detailed script order and input-output map are in `code/README.md`. After
the provider files have been reconstructed, run:

```text
python3 code/30_verify_confirmatory_bundle.py
python3 code/37_verify_dense_bundle.py
python3 code/51_verify_extension_bundle.py
```

The three verifiers check dimensions, simulation sizes, decomposition
identities, retained estimates, extension records, manuscript structure and
the final PDF page counts. They require the full retained simulations rather
than smoke-test output.

To reproduce the added spatial-support, climatology and station checks before
running the verifiers:

```text
python3 code/55_revision_sensitivity_analyses.py
```

## Primary analysis and map display

The primary analysis uses 121 nodes. Script
`31_spatial_field_decomposition.R` calculates exact node allocations for all
five graph scales. Script `38_jrssc_main_figures.R` calls
`53_primary_smooth_spatial_surfaces.R`, which fits low-rank thin-plate REML
surfaces to those values on a 0.1-degree land-masked display grid. The raster
is used only for display; all estimation, uncertainty calculations, tests and
node-sum identities use the 121 observed nodes. The surface and build audit
files record the spline settings, inputs, source-node counts and exact-sum
checks.

## Historical and station extensions

The extension protocol was fixed on 2026-08-09 at 13:50 China Standard Time.
Its SHA-256 is
`da36fc813ffac2dc5e6bab0b316888e598d1e3d02dc00ff70f25fccf77e94251`.
The protocol file is included byte for byte, and the historical, station and
elevation audit files record the same hash.

Scripts 42-44 reconstruct and analyse 121 ERA5-Land points for 41 summers in
1950-1990. The historical energy decomposition estimates its monthly
climatological field within those 41 summers; this decomposition is an
exploratory structural analysis. The dated protocol also specifies the
historical 99,999-draw product shift, so that calculation is a pre-access
extension diagnostic rather than part of the later exploratory global test.
Scripts 44, 45 and 48 produce the continuous-intensity, calendar-matching,
ratio-stress and elevation-basis analyses. Script 46 writes all 108
repeated-summer cells and eight coupled label-field simulation cells to the
portable supplementary table. Script 54 runs the targeted cross-record
dependence size stress test.

Scripts 47, 49 and 50 implement the NOAA comparison for the enumerated set
`{1992+4k: k=0,...,7} union {2023,2025}`. The official
station-history snapshot did not reach the initially specified administrative
end date. The retained, outcome-independent rule uses the latest available
in-scope date in that snapshot, 2025-08-24. The 150-km filter is the minimum
equirectangular distance to the primary 121-site support; it reduces 228
metadata-eligible stations to 208 before the deterministic 30-station
selection. This metadata amendment and the station selection precede retrieval
of station WBT observations, matched ERA5-Land series and graph contrasts.
ERA5-Land supplies finite land-grid
series for 28 stations. The provenance records the two unavailable island
locations, Shengsi (`58472099999`) and Taidong (`59562099999`), and the final
175,172 matched hours. Selection and availability checks do not use graph
effect signs or magnitudes.

## Fixed simulation seeds

- Main day-level graph and variogram simulation: `20260872`
- Repeated-summer simulation: `20260804`
- Coupled label-field simulation: `20260874`
- Irregular 1991-2025 calendar simulation: `20260875`
- Station date-block bootstrap: `20260806`
- Global product cyclic shift: `20260809`
- Historical 1950-1990 product cyclic shift: `20260810`
- Paired heavy-tail ratio experiment: `20260811`
- Cross-record dependence stress test: `20260812`

The retained simulations contain 25,000 main day-level data sets, 8,000
coupled label-field data sets, 8,000 targeted cross-record stress replications,
60,000 irregular-calendar data sets and 1.08 million repeated-summer
replications across 108 cells.

## Runtime recorded for the final analysis

- R 4.3.2
- data.table 1.17.8
- ggplot2 3.5.2
- ggrepel 0.9.6
- patchwork 1.3.1
- maps 3.4.3
- mgcv 1.9-1
- sp 2.2-0
- jsonlite 2.0.0
- ragg 1.4.0
- Python 3.12.1
- NumPy 1.26.4
- pandas 2.2.2
- geopandas 1.0.1
- pyproj 3.6.1
- SciPy 1.14.0
- xarray 2026.7.0
- h5netcdf 1.8.1
- pypdf 6.14.2

The manuscript PDFs were compiled with pdfTeX and LaTeXmk from TeX Live 2024.
