# Extension analysis protocol

Frozen on 2026-08-09 at 13:50 CST, before inspection of any newly acquired
1950-1990 ERA5-Land values or non-development-year NOAA effect estimates.
This protocol adds secondary analyses. It does not redefine the frozen
1991-2025 five-scale estimand.

## Historical ERA5-Land extension

- Period: June-August 1950-1990, all 41 summers.
- Spatial support: the existing 121-site manifest. No site is added, removed,
  or moved.
- Fields: the UTC daily field at the hour of maximum 121-site regional mean,
  constructed from 2-m temperature, 2-m dew-point temperature, and surface
  pressure with the same Bolton wet-bulb implementation and quality rules as
  the frozen analysis.
- Events: type-7 within-month-year quartiles; upper quartile is high, middle
  half is middle, and the lower quartile is unused.
- Graphs: the five fixed physical bandwidths from the primary support,
  125.799765, 251.599530, 503.199060, 1006.398120, and 2012.796241 km.
- Estimator: record-specific high-to-middle ratios, then equal months, equal
  bandwidths, and equal summers.
- Outputs: overall effect, five-scale profile, summer effects, negative-summer
  counts, leave-one-summer-out range, and latitude and planar Laplacian-energy
  decompositions.
- Finite-record diagnostic: one joint product cyclic-shift calculation with
  99,999 draws and seed 20260810. All five dispersion columns move together
  within each month-year record.
- Period check: report 1950-1978 and 1979-1990 separately because ECMWF
  documents a different forcing segment before 1979. This split is based on
  data production, not the observed graph contrast.

The extension will be described as a post-analysis temporal evaluation. The
official ERA5-Land record starts in 1950, so no pre-1950 claim will be made.

## NOAA station extension

- Evaluation years: 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023,
  and 2025. The two development years, 2015 and 2022, are excluded.
- Candidate stations must lie inside the fixed rectangle, be within 150 km of
  an analysis site, and have an ISD history spanning the first through last
  selected summer.
- A deterministic maximin rule selects 30 stations without using WBT effects.
- A station-year enters hourly validation if it has at least 400 exact-hour
  JJA observations with temperature, dew point, and usable station pressure.
  We report the retained station count and hours for every year.
- Measurement summaries are bias, mean absolute error, root mean squared
  error, and within-station centred agreement.
- Scientific-effect summaries use the frozen ERA5-Land peak times and event
  labels. The prespecified station-supported scales are 503, 1006, and 2013
  km; 126 and 252 km remain descriptive.
- Year-specific station subsets are allowed when documented. No station or
  year may be retained or removed according to the sign or magnitude of its
  high-middle effect.

## Structural and sensitivity analyses

- Elevation basis: compare centred latitude, centred latitude-longitude, and
  centred latitude-longitude-elevation column spaces. Elevation comes from an
  official invariant ERA5/ERA5-Land geopotential field at the 121 fixed sites.
  Only total structured and residual energies are reported; no order-specific
  elevation contribution is assigned.
- Date matching: for window widths 3 and 5 days, each high day is compared
  with all middle days in the same month-year whose calendar-day distance is
  within the window. The high-day-specific relative ratios are averaged over
  available high days, then over months, scales, and summers. Coverage and
  records with no eligible match are reported.
- Continuous intensity: within each month-year and bandwidth, regress
  `log(Q_h)` on centred regional-mean WBT. Average slopes equally over months,
  bandwidths, and summers; report scale-specific slopes and negative-summer
  counts. This is descriptive and does not replace the quartile estimand.
- Ratio stress test: reuse the existing heavy-tailed day-level data-generating
  mechanism under its null and -7 percent alternative. Apply raw ratio, log
  ratio, and bounded symmetric contrast to the same simulated fields and the
  same product shifts. Use 2,000 replications per cell, 999 shifts per
  replication, and seed 20260811. Report bias, RMSE, null rejection, and power.
- Existing annual stress test: reproduce all 108 retained cells exactly and
  provide the complete table and data-generating specification in the
  supplement.

## Reporting rules

All new results are labelled post-analysis or temporally separated according
to their actual timing. Existing 1991-2025 labels, bandwidths, and primary
effect measure remain unchanged. Every script records its seed, input hashes,
sample counts, and output paths in a machine-readable audit file.
