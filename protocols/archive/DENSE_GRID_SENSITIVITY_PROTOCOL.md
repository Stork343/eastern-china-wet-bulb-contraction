# Nested-grid spatial-resolution sensitivity protocol

Frozen on 2026-08-03 after completion of the primary 121-site analysis and
before calculation of any 465-site dispersion result. This is a secondary
resolution analysis, not a new confirmatory test. The primary sample, estimand,
and reported inference remain unchanged.

## Purpose and spatial sample

The analysis asks whether the primary spatial-dispersion result is an artifact
of using a visibly sparse set of locations. The domain remains 105--125 degrees
E and 20--42 degrees N. A candidate lattice uses 0.8-degree longitude spacing
and 0.9-degree latitude spacing, with every original site inserted exactly.
The same ERA5-Land reference mask used for the primary design removes ocean
locations. Script `32_define_dense_spatial_grid.R` freezes 465 land sites:
all 121 primary sites and 344 added sites. The original coordinates are never
moved, averaged, or replaced.

The added sites use the same June--August 1991--2025 period and the same
pressure-aware wet-bulb-temperature calculation as the primary analysis. Only
2-m temperature, 2-m dew-point temperature, and surface pressure are acquired
for new sites because the resolution analysis does not repeat the descriptive
mechanism regressions.

## Analysis configurations

The five bandwidths are the exact absolute bandwidths recorded by the primary
analysis, approximately 126, 252, 503, 1,006, and 2,013 km. They are not
re-estimated from the denser sample. Three configurations are calculated:

1. `primary_121_reproduction` repeats the frozen primary computation as a
   numerical quality-control check.
2. `dense_465_fixed_labels` uses the 465-site fields at the peak hours and the
   extreme/moderate labels selected from the original 121-site regional means.
   This isolates the effect of added spatial support.
3. `dense_465_recomputed` selects peak hours and outcome-blind quartile labels
   again from the 465-site regional means. This measures the combined effect of
   added support and the resulting changes in event definition.

The year remains the replication unit. Student, lag-2 Newey--West, and sign
summaries are retained for comparability, but all dense-grid p-values and
intervals are secondary sensitivity results. They do not replace or add to the
frozen three-test confirmation rule.

## Exact node quantities and displayed surfaces

Graph dispersion and its nodewise decomposition are evaluated at the 465
observed locations. Their algebraic identity must hold to numerical tolerance.
No interpolation enters the estimand, graph weights, labels, contrasts, or
uncertainty calculation.

For visualization only, the mean nodewise extreme-minus-moderate contribution
is fitted with a projected-coordinate thin-plate regression spline using REML,
`k = 80`, and shrinkage selection. Predictions use a 0.1-degree display grid.
Pixels are retained only on land and within 140 km of an observed site. Small
site markers are drawn over the surface so that readers can see both the data
support and the descriptive interpolation. The manuscript must call these
surfaces smoothed displays, not estimated continuous spatial effects.

## Quality gates

The dense-grid result may be reported only if all of the following checks pass:

1. the manifest contains 465 unique ERA5-Land coordinates and nests the 121
   primary coordinates exactly;
2. all 344 new point archives cover the requested hourly period and their
   stored coordinates match the manifest;
3. every year contains 92 complete daily fields at all 465 sites;
4. the 121-site reproduction differs from the frozen result by no more than
   `1e-9` at every matched year and scale;
5. the exact nodewise decomposition identity has error no greater than
   `1e-10`; and
6. the fixed-label and recomputed-label analyses, their audits, and both smooth
   surface figures are present in the completion bundle.

After the 465-site results are calculated, changes to these choices are limited
to documented software corrections. Any substantive alternative must be
reported separately as exploratory.
