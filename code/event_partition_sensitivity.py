"""Post-decision event-partition sensitivity; primary definitions stay fixed."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from extended_analyses import load_primary_fields, projected_distance, edge_field_metrics, H_FACTORS, t_summary

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/event_partition'
OUT.mkdir(parents=True, exist_ok=True)
sites, days, matrix = load_primary_fields()
distance = projected_distance(sites)
bandwidths = np.median(distance[np.triu_indices(len(sites), 1)]) * H_FACTORS
q = edge_field_metrics(matrix, distance, bandwidths)
variance = np.var(matrix, axis=0)
metrics = np.column_stack([q, variance])
specs = [('Quartiles', .25, .75), ('Tertiles', 1/3, 2/3), ('20--80', .2, .8), ('30--70', .3, .7), ('25--70', .25, .7), ('25--80', .25, .8)]
rows = []
for name, lower, upper in specs:
    for record, indices in days.groupby('record_id').indices.items():
        block = days.iloc[indices]
        year = int(block.year.iloc[0])
        if year in (2015, 2022):
            continue
        mean = block.regional_mean_wbt.to_numpy()
        lo, hi = np.quantile(mean, [lower, upper], method='linear')
        high, middle = mean >= hi, (mean >= lo) & (mean < hi)
        assert high.any() and middle.any()
        a, b = metrics[indices][high].mean(axis=0), metrics[indices][middle].mean(axis=0)
        for k in range(6):
            rows.append(dict(partition=name, lower=lower, upper=upper, record=record, year=year, scale=k+1, high_n=int(high.sum()), middle_n=int(middle.sum()), ratio=a[k]/b[k]-1, rms_absolute=np.sqrt(2*a[k])-np.sqrt(2*b[k])))
records = pd.DataFrame(rows)
yearly = records.groupby(['partition', 'year', 'scale'], as_index=False)[['ratio', 'rms_absolute']].mean()
summary = []
for name, _, _ in specs:
    sub = yearly[yearly.partition == name]
    for scale in range(0, 7):
        values = sub[sub.scale <= 5].groupby('year').ratio.mean() if scale == 0 else sub[sub.scale == scale].ratio
        summary.append(dict(partition=name, scale=scale, **t_summary(values.to_numpy())))
summary = pd.DataFrame(summary)
expected = pd.read_csv(ROOT / 'results/confirmatory_primary_results.csv')
print(expected[['day_definition','estimate']].to_string(index=False))
primary = summary[(summary.partition == 'Quartiles') & (summary.scale == 0)].estimate.iloc[0]
error = abs(primary - expected[expected.day_definition == 'utc'].estimate.iloc[0])
assert error < 1e-10, error
records.to_csv(OUT/'record_contrasts.csv',index=False)
yearly.to_csv(OUT/'year_contrasts.csv',index=False)
summary.to_csv(OUT/'summary.csv',index=False)
(OUT/'audit.json').write_text(json.dumps(dict(status='post-decision exploratory sensitivity', date='2026-09-07', primary_reproduction_error=error, records_per_partition=99, years=33, quantile='type 7 / linear', intervals='descriptive t-reference; no serial-dependence coverage guarantee'),indent=2))
print(summary[summary.scale.isin([0,1,5,6])].to_string(index=False))
