# September 2026 application validation

This finite extension uses the original 121-site daily peak and event labels,
465-site `primary_grid_peak` fields, all 33 evaluation summers and a
77-threshold development-defined grid. The methods and stopping rule are in
`output_revision_application/revision_protocol.md` and `config.json`.

For a fresh recomputation, preserve the supplied output folder as
`output_revision_application_reference` and let step 1 create a new
`output_revision_application`. The versioned protocol is copied from this
code directory before outcome calculations.

Run from the project root with Python 3.12 (NumPy, pandas, SciPy, scikit-learn,
pyarrow, geopandas, shapely, pyproj and matplotlib). R with spData supplies the
existing Natural Earth land boundary. Input acquisition is unchanged.

1. `python3 code/revision_202609/00_audit.py` (fresh output directory only;
   refuses to overwrite a frozen config).
2. `python3 code/revision_202609/01_footprints.py`
3. `python3 code/revision_202609/02_regime.py`
4. `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/revision_202609/03_compare.py`
5. `python3 code/revision_202609/04_summarise.py`
6. `python3 code/revision_202609/05_verify.py` before editing the manuscript.
7. `python3 code/revision_202609/06_figures_tables.py`

Existing source data, original code, protocols and results are protected by
SHA-256 comparisons. The manuscript changes follow checkpoint3. The figure
and table generator rebuilds numeric tables and figures; the narrative
LaTeX snippets in output_revision_application/generated are author-facing
manuscript revisions rather than numerical output.

The regression is weighted ridge least squares for bounded proportions, with
quadratic splines and [0,1] clipping. Penalties use nested blocked validation.
It is an omitted-years comparison within one reanalysis system, not a weather
forecast or an independent measurement validation. Interpolation preserves
full primary-field information and uses no dense temperatures as inputs.

Area targets approximate the retained land at 0.8° × 0.9° sampling resolution.
Coverage curves use strict exceedance. Only fixed grid neighbours determine
connected area. The additional 344-site check is for area, not connectivity.

## Formatting the complete simulation tables

After regenerating the retained simulation tables, create independent tables
for each sample size with:

```text
python3 code/revision_202609/09_paginate_tables.py --input results/supp_complete_simulation_tables.tex --output manuscript/generated/supp_complete_simulation_tables.tex
```

These paths refer to the public archive. In the journal workspace, use
`JRSSC/manuscript/generated/supp_complete_simulation_tables.tex` as the output.
The formatter preserves every numerical row and does not run simulations.
