"""Paired annual held-out MAE, full threshold curves and fixed sensitivities."""
import numpy as np
import pandas as pd
from common import OUT,config,save_json,verify_originals,sha256
conf=config();f=pd.read_parquet(OUT/'heldout_predictions.parquet')
assert not f.duplicated(['method','analysis_date','threshold_c']).any()
assert len(f)==5*3036*len(conf['thresholds_c'])
rows=[]
for subset in ['all','high']:
 d=f if subset=='all' else f[f.regime=='high']
 for endpoint in ['A','C','A_new344']:
  t=d[['year','month','method','threshold_c']].copy()
  t['mae']=abs(d[endpoint+'_predicted']-d[endpoint+'_observed'])
  t['bias']=d[endpoint+'_predicted']-d[endpoint+'_observed']
  a=t.groupby(['year','month','method','threshold_c'],as_index=False)[['mae','bias']].mean().groupby(['year','method','threshold_c'],as_index=False)[['mae','bias']].mean()
  a['outcome']=endpoint;a['subset']=subset;rows.append(a)
curve=pd.concat(rows,ignore_index=True);curve.to_csv(OUT/'annual_threshold_comparison.csv',index=False)
annual=curve.groupby(['year','method','outcome','subset'],as_index=False)[['mae','bias']].mean();annual['outer_fold']=(annual.year-1991)//5+1
annual.to_csv(OUT/'annual_comparison.csv',index=False)
summary=annual.groupby(['method','outcome','subset'],as_index=False)[['mae','bias']].mean();summary['mae_pp']=100*summary.mae;summary.to_csv(OUT/'comparison_summary.csv',index=False)
curve.groupby(['method','outcome','subset','threshold_c'],as_index=False)[['mae','bias']].mean().to_csv(OUT/'threshold_comparison.csv',index=False)
boot=np.load(OUT/'calendar_bootstrap_weights.npy')
pairrows=[];pairstats=[];paircurves=[]
for (endpoint,subset),a in annual.groupby(['outcome','subset']):
 p=a.pivot(index='year',columns='method',values='mae').sort_index()
 for ref in ['simple','variogram']:
  diff=p[ref]-p.graph
  for year,value in diff.items():pairrows.append(dict(year=year,outcome=endpoint,subset=subset,reference=ref,improvement_pp=100*value))
  draws=boot@diff.to_numpy();lo,hi=np.quantile(draws,[.025,.975])
  pairstats.append(dict(outcome=endpoint,subset=subset,reference=ref,improvement_pp=float(100*diff.mean()),lower_pp=100*lo,upper_pp=100*hi,improved_years=int(sum(diff>0)),worse_years=int(sum(diff<0)),relative_reduction=float(diff.mean()/p[ref].mean())))
  curves=curve[(curve.outcome==endpoint)&(curve.subset==subset)].pivot(index=['year','threshold_c'],columns='method',values='mae')
  vv=(curves[ref]-curves.graph).unstack('threshold_c').sort_index()
  lo,hi=np.quantile(boot@vv.to_numpy(),[.025,.975],axis=0)
  for k,u in enumerate(vv.columns):paircurves.append(dict(outcome=endpoint,subset=subset,reference=ref,threshold_c=u,improvement_pp=100*vv.iloc[:,k].mean(),lower_pp=100*lo[k],upper_pp=100*hi[k]))
pd.DataFrame(pairrows).to_csv(OUT/'annual_paired_improvement.csv',index=False)
pd.DataFrame(pairstats).to_csv(OUT/'paired_improvement_summary.csv',index=False)
pd.DataFrame(paircurves).to_csv(OUT/'threshold_paired_improvement.csv',index=False)
# Preserve both alternating coarser threshold grids, with no refitting.
u=np.array(conf['thresholds_c']);coarse=[]
for offset in [0,1]:
 a=curve[curve.threshold_c.isin(u[offset::2])].groupby(['year','method','outcome','subset'],as_index=False).mae.mean();a['offset']=offset;coarse.append(a)
pd.concat(coarse).to_csv(OUT/'threshold_grid_sensitivity.csv',index=False)
fold=pd.read_csv(OUT/'fold_manifest.csv')
for outer in range(1,8):
 o=fold[(fold.outer_fold==outer)&(fold.inner_fold==0)];test=o[o.role=='test'].year.to_numpy();tr=o[o.role=='train'].year.to_numpy()
 assert set(test).isdisjoint(tr) and np.min(abs(test[:,None]-tr[None,:]))>=2
 for inner,ii in fold[(fold.outer_fold==outer)&(fold.inner_fold>0)].groupby('inner_fold'):
  va=ii[ii.role=='validation'].year.to_numpy();itr=ii[ii.role=='train'].year.to_numpy()
  assert set(itr).issubset(set(tr)) and set(va).issubset(set(tr)) and np.min(abs(va[:,None]-itr[None,:]))>=2
assert len(f[f.method=='graph'].analysis_date.unique())==3036
# Required checkpoint reached before any manuscript edits.
save_json('checkpoint3.json',dict(status='passed',heldout_rows=len(f),methods=5,outcomes=3,evaluation_days=3036,thresholds=len(u),fold_separation_verified=True,config_unchanged=sha256(OUT/'config.json')==__import__('json').loads((OUT/'checkpoint0.json').read_text())['config_sha256'],protocol_unchanged=sha256(OUT/'revision_protocol.md')==__import__('json').loads((OUT/'checkpoint0.json').read_text())['protocol_sha256'],original_files=verify_originals()))
print(summary.pivot(index='method',columns=['outcome','subset'],values='mae_pp').to_string())
print(pd.DataFrame(pairstats).to_string(index=False))
