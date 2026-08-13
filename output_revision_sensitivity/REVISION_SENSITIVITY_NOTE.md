# Revision sensitivity analysis note

These are post-review sensitivity analyses; they do not replace the protocol-defined primary finite-record summary. All reported intervals are t-scaled between-summer variability intervals, not confidence intervals under serial dependence.

## Spatial support and distance

The frozen equirectangular/equal-site calculation reproduced the stored result to 5.00e-16. Changing only distance to WGS84 gave -6.97%, versus -7.28% for the original calculation. WGS84/cosine-latitude weighting gave -7.80% with frozen labels and -7.37% after relabelling.

The non-rectangular land target uses user-supplied Natural Earth-compatible file: the country feature named China, excluding the separately stored Taiwan feature, intersected with the downloaded dense lattice and its frozen ERA5-Land-valid mask. It retained 445 sites and gave -6.35%. All one-edge sensitivity analyses move the edge inward; outward movement is not estimable without new ERA5-Land downloads.

## Leave-one-summer-out climatology

The exact Eq. 11 leave-one-summer-out decomposition gave -1.27 percentage points for anomaly energy and -6.01 percentage points for the climatology--anomaly cross term on the five-bandwidth average. This is an algebraic structural interpretation, not a causal mechanism attribution.

## NOAA station comparison

Frozen ERA labels with dynamic support gave -17.37% on the broad-scale average. A year-specific common-station support gave -18.08%, and station-defined peak hours and event labels with at least five days per regime gave -15.99%.

Among 680 frozen high/middle event fields, 345 had at least 10 stations. The high-minus-middle mean station-count difference was 0.16 stations. The fixed-support analysis uses a separate station intersection within each summer; the all-summer intersection contains only 2 stations and cannot support a regional graph.

Station pressure was computed as `p_station_Pa = SLP_hPa * 100 * max(1 - 2.25577e-5 * elevation_m, 0.1)^5.2559`. Station elevations ranged from 3 to 2063 m. Treating sea-level pressure directly as surface pressure changed WBT by 0.072 degrees C on average in absolute value and at most 1.346 degrees C.
