JRSS C source package
=====================

Main source: main.tex
Supplementary source: supplement_theory.tex
Official OUP class: oup-authoring-template.cls (version 1.5, 2026-07-14)
OUP bibliography styles: oup-abbrvnat.bst and oup-plain.bst
Figures: figures/*.pdf (seven main-text and two supplementary vector PDFs)
Generated table: generated/supp_complete_simulation_tables.tex

Compile the main manuscript from the package root with:

    latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

Compile the supplement with:

    latexmk -pdf -interaction=nonstopmode -halt-on-error supplement_theory.tex

Both documents use the official OUP modern, medium, single-column template
with author-year citations. The bibliography is embedded in main.tex, and all
paths resolve inside this archive. The generated file contains all 108 repeated-summer cells and
the eight coupled label-field simulation cells. The supplementary source also
reports the four-cell cross-record dependence size stress test; its code and
machine-readable results are in the separate scientific reproducibility bundle.

Figures 3 and 4 display smooth surfaces fitted to values from the primary
121-node analysis. All 121 source nodes are shown. Raster cells are graphical
interpolations and do not enter estimation, uncertainty calculations, tests or
node-sum identities.

The compiled supplementary PDF is below the JRSS C limit of 2 MB per
supplementary-material file. The larger scientific reproducibility archive is
frozen at
https://github.com/Stork343/eastern-china-wet-bulb-contraction/releases/tag/jrssc-submission-v1-2026-08-13.
