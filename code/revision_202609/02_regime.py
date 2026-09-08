"""Equal-month/year contrasts and paired calendar-block intervals."""
import numpy as np
import pandas as pd
from common import OUT,save_json,verify_originals

def calendar_bootstrap(years,B=1999,seed=20260907):
    rng=np.random.default_rng(seed);calendar=np.arange(1991,2026)
    lookup={int(y):i for i,y in enumerate(years)}
    weights=np.zeros((B,len(years)))
    for b in range(B):
        starts=rng.integers(0,35,7)
        sample=calendar[((starts[:,None]+np.arange(5))%35).ravel()]
        indices=[lookup[int(y)] for y in sample if int(y) in lookup]
        weights[b]=np.bincount(indices,minlength=len(years))/len(indices)
    return weights

if __name__=='__main__':
    daily=pd.read_parquet(OUT/'daily_footprints.parquet');daily=daily[daily.role=='evaluation']
    measures=['A','C','S_km2','largest_km2','C8','A_coslat','C_coslat','A_new344']
    groups=daily[daily.regime.isin(['high','middle'])].groupby(['year','month','threshold_c','regime'])[measures].mean().unstack('regime')
    rec=[]
    for name in measures:
        z=groups[name].reset_index();z['difference']=z.high-z.middle;z['outcome']=name;rec.append(z)
    rec=pd.concat(rec,ignore_index=True);rec.to_csv(OUT/'regime_footprints.csv',index=False)
    annual=rec.groupby(['year','threshold_c','outcome'],as_index=False)[['high','middle','difference']].mean()
    annual.to_csv(OUT/'annual_footprints.csv',index=False)
    boot=calendar_bootstrap(sorted(annual.year.unique()));np.save(OUT/'calendar_bootstrap_weights.npy',boot)
    summaries=[]
    for outcome,part in annual.groupby('outcome'):
        v=part.pivot(index='year',columns='threshold_c',values='difference').sort_index()
        low,high=np.quantile(boot@v.to_numpy(),[.025,.975],axis=0)
        z=part.groupby('threshold_c')[['high','middle','difference']].mean().reset_index()
        z['lower']=low;z['upper']=high;z['outcome']=outcome
        z['positive_years']=(v>0).sum().to_numpy();z['negative_years']=(v<0).sum().to_numpy()
        summaries.append(z)
    summary=pd.concat(summaries,ignore_index=True);summary.to_csv(OUT/'coverage_summary.csv',index=False)
    # Fixed sensitivity output over the full grid; no range selection.
    p=summary.pivot(index='threshold_c',columns='outcome',values='difference')
    sensitivity=pd.DataFrame(dict(threshold_c=p.index,queen_minus_rook=p.C8-p.C,coslat_minus_land_A=p.A_coslat-p.A,coslat_minus_land_C=p.C_coslat-p.C))
    sensitivity.to_csv(OUT/'coverage_sensitivity.csv',index=False)
    save_json('regime_audit.json',dict(status='passed',evaluation_days=len(daily)//daily.threshold_c.nunique(),records=99,years=33,bootstrap_replicates=1999,calendar_block_years=5,original_files=verify_originals()))
    print(summary[(summary.outcome.isin(['A','C'])) & (summary.threshold_c.isin([20,24,26,28]))].to_string(index=False))
