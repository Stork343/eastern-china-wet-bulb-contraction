# Confirmatory analysis protocol

Version 2, frozen on 2026-08-02 before acquisition or inspection of the
multi-year ERA5-Land point panel. Version 2 replaces the single-test rule in
version 1 after outcome-free stress simulations showed that no candidate
year-level test was reliable enough on its own. No confirmatory ERA5-Land
outcome was available when this change was made.

## Scope and analysis roles

- Study window: June-August 1991-2025.
- Discovery years: 2015 and 2022. These years were used to develop and correct
  the method and cannot contribute to the primary confirmatory test.
- Confirmatory years: the other 33 years in 1991-2025.
- The primary analysis uses UTC days. Complete UTC+8 calendar days form a
  prespecified sensitivity analysis.
- All 35 years will be shown in descriptive figures, but confirmatory estimates
  and p-values will be visibly separated from discovery results.

## Spatial sample

The domain is 105-125 degrees E and 20-42 degrees N on the 0.1-degree
ERA5-Land grid. Starting at the north-west corner, the candidate lattice takes
every 16th longitude cell and every 18th latitude cell. Candidate locations
with missing ERA5-Land 2-m temperature in the cached 2015-06 reference field
are treated as ocean and removed. Script `20_define_spatial_grid.R` reproduces
the resulting 121 locations exactly and freezes their coordinates in
`data/grid/eastern_china_121_sites.csv`.

## Daily fields and quality control

For every day, the primary field is the simultaneous 121-site vector at the
hour of maximum regional mean pressure-aware WBT. Ties are resolved by the
earliest hour. Required quality gates are:

1. 121 finite values at every retained hour;
2. 2-m temperature between 220 and 340 K;
3. surface pressure between 45,000 and 110,000 Pa;
4. dew point clipped to air temperature only when interpolation produces a
   positive exceedance, with the count reported by year;
5. exactly 92 daily fields per summer and day definition; and
6. no change in site identifiers or coordinates across years.

The UTC+8 sensitivity uses a 24-hour acquisition buffer so that 1 June and
31 August are complete local calendar days.

## Outcome-blind classification

Classification is performed separately within each month-year record. A day
is extreme when regional mean WBT is at or above its record-specific 75th
percentile. Moderate days lie between the 25th and 75th percentiles. Lower
quartile days are excluded. No graph dispersion, spatial variance, or
mechanism variable enters the labels.

## Fixed spatial outcomes

Graph weights are Gaussian functions of projected distance. Five bandwidths
are fixed at 0.125, 0.25, 0.5, 1, and 2 times the median pairwise distance of
the frozen grid. They correspond to approximately 126, 252, 503, 1,006, and
2,013 km. For each scale, the outcome is normalized graph dispersion

\[
Q_h(y)=\frac{y^T L_h y}{2\sum_{i<j}w_{ij,h}}.
\]

The scale-specific contrast is the extreme-to-moderate mean ratio minus one.
The primary profile statistic is the unweighted mean of the five relative
contrasts. Spatial variance is a prespecified complete-graph benchmark.

## Primary estimand and inference

For every confirmatory year, each scale-specific relative contrast is first
calculated within June, July, and August and then averaged equally across the
three months. The five scale effects are averaged to obtain one profile effect
per year. The primary estimand is the mean of these 33 yearly profile effects.
A related recurrence estimand is the fraction of confirmatory years with a
negative profile effect.

The confirmatory claim requires all three prespecified one-sided results to be
at or below 0.025: (i) the ordinary Student test of the mean, (ii) a Newey-West
test of the mean with lag 2, and (iii) the exact binomial sign test of whether
negative yearly effects recur more often than chance. The first two address
effect magnitude under different year-dependence assumptions; the third
addresses direction recurrence rather than the mean. Two-sided 95% ordinary
and Newey-West intervals, the number of negative years, and leave-one-year-out
means will all be reported. The composite rule deliberately requires agreement
about both mean magnitude and typical direction; failure of any component is
reported as a failure to confirm the claim.

In 10,000 outcome-free simulations with 33 yearly effects, the original
Student test rejected 9.1% of skewed mean-zero samples and 11.0% of AR(0.3)
samples at nominal 0.05. The three-test 0.025 consistency rule rejected 0.5%
to 4.6% across Gaussian, skewed, heavy-tailed, single-outlier, and AR(0.3)
null scenarios. These simulations motivated the revision; they are not fitted
to the observed multi-year result. The year, rather than the day, remains the
replication unit.

Scale-specific one-sided Student, Newey-West, and sign-test values are each
adjusted by Holm's procedure across the five fixed bandwidths. A scale is
called confirmed only if all three adjusted values are at or below 0.025.
Their confidence intervals and unadjusted values will also be reported.

## Secondary analyses

- Spatial variance and a within-record continuous mean-WBT slope.
- UTC+8 day definition.
- Extreme thresholds at the 70th and 80th percentiles.
- Stull WBT and sitewise daily maxima.
- Early (1991-2007) versus late (2008-2025) effect estimates, treated as a
  heterogeneity analysis rather than evidence of a trend.
- Descriptive associations with 10-m wind, surface solar radiation, top-layer
  soil moisture, and separately acquired ERA5 pressure-level circulation
  fields.

## Reporting restrictions

- Results from 2015 and 2022 will not be called validation results.
- Scale choices and the negative direction will be identified as developed in
  the discovery analysis.
- Mechanism covariates may support consistency with a physical explanation but
  will not establish causality.
- Any deviation from this protocol will be dated, justified, and reported in a
  separate deviations section before the affected result is calculated.

## Recorded implementation deviations (2026-08-02)

1. The official CDS ERA5-Land time-series API replaced direct ARCO Zarr access.
   Direct Zarr metadata and chunk loading exceeded practical memory and elapsed
   time on the analysis machine. Both routes expose ERA5-Land, and the API
   requests used the frozen coordinates, dates, and variables without spatial
   interpolation. This change was made before confirmatory outcomes were
   analyzed. Raw point archives, request metadata, and SHA-256 hashes are
   retained.
2. NetCDF packing yielded 2,092 hourly dew-point values slightly above air
   temperature among 9,452,520 time-site records (0.022%). The protocol
   anticipated positive exceedances after interpolation, but the point API
   required no interpolation. The same physical consistency rule was applied:
   dew point was replaced by air temperature at those records, with yearly
   counts written to the audit table. No outcome-dependent threshold or record
   exclusion was introduced.
3. ERA5 pressure-level circulation fields were not acquired in this run.
   Mechanism summaries are therefore restricted to the prespecified ERA5-Land
   wind, surface solar radiation, and top-layer soil-water variables and remain
   descriptive.
