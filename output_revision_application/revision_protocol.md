# Simultaneous humid-heat coverage: analysis protocol

This application extension follows the September 2026 manuscript review. All earlier primary, historical, station and simulation results remain unchanged. The protocol and `config.json` are fixed before calculating evaluation-period coverage contrasts or fitting comparison models.

## Question and outcomes

At the original daily regional peak time, how much land simultaneously exceeds wet-bulb temperature u, and how much belongs to the largest connected exceedance region? How well do statistics from 121 sites describe these quantities on the 465-site support?

Use `analysis_definition == primary_grid_peak` exclusively and the original 121-site UTC peak, type-7 high/middle/low labels, 33 evaluation summers and five physical bandwidths. Development summers 2015 and 2022 select the temperature range and are excluded from reported evaluation contrasts and model fitting. An additional analysis remains exploratory despite being specified before its own calculation.

At each threshold use strict exceedance `wbt > u`. A is weighted exceedance proportion; S is A times represented land area in km². C is largest connected exceedance component area divided by total represented land area. Zero exceedance gives A=C=0. Four-neighbour grid connectivity is primary and eight-neighbour connectivity is the sole connectivity sensitivity. Gaussian graph edges are never used for connectivity. All sites in a field share one observed time.

## Spatial support

Representative cells extend half the 0.8° longitude and 0.9° latitude spacing around each of the 465 retained sites, clipped to the original rectangle 105–125°E, 20–42°N and a Natural Earth land polygon. Use the land union across all countries within that rectangle, matching the geographical support of the existing analysis, rather than selecting a new China-only subset. Compute land area on WGS84. Omitted ERA5-Land-invalid cells remain outside the denominator. Report total represented land and its fraction of the rectangle's mapped land. Coastline detail and connected areas remain approximations at the sampling resolution. A partially coastal cell has a single observed value. Grid neighbours describe cell-scale connectivity, not subcell coastal corridors. Zero-area cells cannot bridge components.

The complete 465-point reference includes the original 121 points. A second exceedance endpoint uses only the 344 additional points and renormalises their land areas. C always uses the entire 465-point grid. The 121-point direct exceedance feature uses static representative weights obtained by assigning each dense-cell area to its nearest primary site in the existing projected coordinates; this is a discretised Voronoi approximation on the same area support.

## Temperature thresholds

Pool the synchronous 465-site values from the two development summers only. Compute their 1st and 99.5th sample percentiles, round the lower endpoint downward and the upper endpoint upward to multiples of 0.25°C, and include every 0.25°C threshold between them. Record exact endpoints in config.json. Report all thresholds. These are temperature levels, not clinical warning thresholds. Check integrated errors on the two alternating 0.5°C subgrids without refitting or selecting one.

## Coverage comparison

For every month-year, compute high-minus-middle means of A, S, C and the largest-component area in km². Average June, July and August equally, then the 33 summers equally. Report percentage-point and km² differences. Retain both group means and every annual threshold curve, including zero crossings. Calculate full four/eight-neighbour and land-area/cosine-latitude comparisons once; no factorial search over old model choices.

For descriptive uncertainty, resample five-year blocks of the calendar sequence 1991–2025, retaining the positions of the two excluded development years. Draw 1,999 circular moving-block replicates with seed 20260907, truncate to 35 calendar positions, then drop development-year positions and average the available evaluation summers. Use shared draws across all thresholds and paired methods. Report pointwise percentile bands, with no simultaneous significance claim. Daily fields are not independent replicates.

## Fixed model comparison

Inputs come only from the original 121-site field. Compare these five methods:

1. Regional mean and day-of-summer.
2. The same features plus ordinary spatial variance, signed OLS latitude gradient in °C per degree latitude, and primary-site direct exceedance proportion.
3. The strong simple baseline plus five log(binned semivariance / variance) values. Distance bins are fixed pair-distance quintiles of the original coordinates, with all 7,260 pairs assigned once.
4. The strong simple baseline plus five log(Q_h / variance) values, using original physical bandwidths.
5. Piecewise-linear interpolation of the 121 temperatures to the dense grid in the original projected kilometre coordinates, with nearest-site extrapolation beyond the convex hull, then direct coverage calculation. This is a full-field reference.

For the four summary models use the same weighted ridge least-squares proportion regression. Give each baseline scalar a quadratic B-spline basis with four uniformly spaced knots fitted on training data, including day-of-summer; add graph and variogram shape features linearly. Standardise training features, fit an unpenalised intercept and penalised slopes, and clip predicted proportions to [0,1]. This fixes model complexity and avoids treating weighted proportions as counts of independent Bernoulli sites. There is no monotonicity or C≤A postprocessing for regression predictions; physical identities apply to measured fields. Report any inconsistent predictions as an approximation diagnostic.

If V=0, define all log-shape features as zero by convention and record the count. If V>0 and any graph/bin energy is nonpositive, stop and diagnose rather than adding an arbitrary epsilon.

Fit A, C and additional-344-point A at each threshold. Select one ridge penalty per method, endpoint and outer fold from {0.1,1,10,100,1000}, minimising inner held-out MAE averaged equally over thresholds, months and years. Use identical processing and candidate penalties for all summary methods. Sample weights give equal weight to each month in each training year.

## Temporal validation

Use outer calendar blocks 1991–1995, 1996–2000, …, 2021–2025. Hold out all evaluation summers in one block and exclude adjacent calendar years from training. Omit development summers throughout fitting. Inner folds are the remaining original calendar blocks, restricted to outer training years, again with a one-year buffer. Fit all splines, standardisation and penalties within the relevant training fold. The time split assesses transport across omitted years, not forecasting: training years can occur on both sides of the held-out period. Store outer and inner memberships and selected penalties.

For each method and summer report MAE integrated equally over all thresholds, first averaging within month, then over months. Repeat on high-mean days. Also retain threshold-specific losses and predictions for every held-out day. Primary comparisons are strong simple minus graph and variogram minus graph, in percentage points; retain negative improvements and all years. Interval resampling of annual paired losses conditions on the fitted cross-validation predictions and does not refit the models or imply independent fold estimates.

## Endpoints and stopping rule

Checkpoint 0: data audit passes, weights and adjacency are interpretable, protocol/config hashes recorded. Proceed without using outcome differences to alter the design.

Checkpoint 3: all five methods have predictions for every held-out evaluation day and all thresholds; paired annual errors and temperature-specific changes are available. Interpret the results, then revise the paper once. No new large simulations, model roster expansion, favourable-threshold selection or historical/station reanalysis. Corrections to implementation errors are recorded.

Outputs follow the filenames requested by the authors under `output_revision_application/`. Code is in `code/revision_202609/`.

## Sources informing the design

- JRSS C author instructions: https://academic.oup.com/jrsssc/pages/general-instructions
- Healy et al., spatial extreme-temperature events: https://academic.oup.com/jrsssc/article/74/2/275/7749363
- Spatial extent of heat waves: https://www.nature.com/articles/s43247-025-02661-y
- Roberts et al., cross-validation for dependent data: https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/ecog.02881
