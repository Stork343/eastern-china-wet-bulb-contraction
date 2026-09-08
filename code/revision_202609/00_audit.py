"""Audit fixed synchronous supports and freeze application comparisons before outcomes."""
import json, math, subprocess
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from shapely.ops import unary_union
from pyproj import Geod
from common import ROOT,OUT,CODE,DEVELOPMENT,sha256,save_json
OUT.mkdir(exist_ok=True)
if (OUT/'config.json').exists():raise RuntimeError('Protocol already frozen; preserve supplied outputs and start in a fresh output directory')
# Restore the versioned design before any derived outcomes on a fresh rerun.
import shutil
shutil.copy2(CODE/'revision_protocol.md',OUT/'revision_protocol.md')
protected={}
for prefix in ['results','JASA']:
 for p in sorted((ROOT/prefix).rglob('*')):
  if p.is_file() and p.name!='.DS_Store':protected[str(p.relative_to(ROOT))]=sha256(p)
for p in sorted((ROOT/'code').glob('*')):
 if p.is_file():protected[str(p.relative_to(ROOT))]=sha256(p)
for p in ROOT.glob('*PROTOCOL.md'):protected[p.name]=sha256(p)
save_json('original_hashes.json',protected)
save_json('manuscript_before.json',{str(p.relative_to(ROOT)):sha256(p) for p in (ROOT/'JRSSC/manuscript').glob('*.tex')})
grid=pd.read_csv(ROOT/'data/grid/eastern_china_dense_sites.csv').sort_values('dense_site_id').reset_index(drop=True)
primary=pd.read_csv(ROOT/'data/grid/eastern_china_121_sites.csv').sort_values('site_id')
assert len(grid)==465 and len(primary)==121
key=lambda x:list(zip(np.round(x.lon,6),np.round(x.lat,6)))
lookup={k:i for i,k in enumerate(key(grid))}
pidx=np.array([lookup[k] for k in key(primary)])
assert len(set(pidx))==121
assert np.array_equal(grid.iloc[pidx].original_site_id.to_numpy(),primary.site_id)
grid['row']=np.rint((grid.lat-20.4)/.9).astype(int)
grid['col']=np.rint((grid.lon-105)/.8).astype(int)
assert np.allclose(grid.lat,20.4+.9*grid.row) and np.allclose(grid.lon,105+.8*grid.col)
assert not grid.duplicated(['row','col']).any()
boundary=subprocess.check_output(['Rscript','-e','cat(system.file("shapes/world.gpkg",package="spData"))'],text=True).strip()
world=gpd.read_file(boundary).to_crs(4326)
region=box(105,20,125,42)
land=unary_union(list(world.geometry)).intersection(region)
geod=Geod(ellps='WGS84')
# Preserve holes and polygon orientation for ellipsoidal land areas.
from shapely.geometry.polygon import orient
def area(g):
 if g.is_empty:return 0.
 if g.geom_type=='Polygon':return abs(geod.geometry_area_perimeter(orient(g,sign=1))[0])/1e6
 return sum(area(v) for v in getattr(g,'geoms',[]))
cells=[box(r.lon-.4,r.lat-.45,r.lon+.4,r.lat+.45).intersection(region) for r in grid.itertuples()]
clips=[c.intersection(land) for c in cells]
grid['land_area_km2']=[area(c) for c in clips]
grid['cell_area_km2']=[area(c) for c in cells]
grid['land_fraction']=grid.land_area_km2/grid.cell_area_km2
assert (grid.land_area_km2>=0).all() and (grid.land_area_km2<=grid.cell_area_km2*1.001).all()
grid['weight']=grid.land_area_km2/grid.land_area_km2.sum()
grid['coslat_weight']=np.cos(np.deg2rad(grid.lat));grid.coslat_weight/=grid.coslat_weight.sum()
grid.to_csv(OUT/'cell_weights.csv',index=False)
gpd.GeoDataFrame(grid,geometry=clips,crs=4326).to_file(OUT/'represented_land_cells.gpkg',driver='GPKG')
# Adjacency is sampling-grid adjacency; absent or zero-area nodes cannot bridge components.
rc={(r.row,r.col):i for i,r in enumerate(grid.itertuples()) if r.land_area_km2>0}
edges=[]
for (row,col),i in rc.items():
 for dr,dc in [(0,1),(1,0),(1,1),(1,-1)]:
  j=rc.get((row+dr,col+dc))
  if j is not None:edges.append(dict(i=i,j=j,rook=(dr==0 or dc==0)))
pd.DataFrame(edges).to_csv(OUT/'adjacency.csv',index=False)
manifest=[];pfields=[];dfields=[];audits=[];inputs={}
for year in range(1991,2026):
 pp=ROOT/f'data/era5_confirmatory/daily_fields/era5_land_{year}_jja_daily_fields.csv.gz'
 dp=ROOT/f'data/era5_dense/daily_fields/era5_land_{year}_jja_dense_daily_fields.csv.gz'
 inputs[str(pp.relative_to(ROOT))]=sha256(pp);inputs[str(dp.relative_to(ROOT))]=sha256(dp)
 p=pd.read_csv(pp);d=pd.read_csv(dp)
 p=p[p.day_definition=='utc'].copy();d=d[d.analysis_definition=='primary_grid_peak'].copy()
 dates=pd.date_range(f'{year}-06-01',f'{year}-08-31').strftime('%Y-%m-%d').tolist()
 assert len(p)==92*121 and len(d)==92*465
 for frame,n in [(p,121),(d,465)]:
  assert not frame.duplicated(['analysis_date','site_id']).any()
  assert sorted(frame.analysis_date.unique())==dates
  assert (frame.groupby('analysis_date').site_id.nunique()==n).all()
  assert (frame.groupby('analysis_date').time.nunique()==1).all()
  coords=frame[['site_id','requested_lon','requested_lat']].drop_duplicates().sort_values('site_id')
  expect=primary if n==121 else grid
  assert len(coords)==n and np.allclose(coords[['requested_lon','requested_lat']],expect[['lon','lat']],atol=1e-8)
 a=p.pivot(index='analysis_date',columns='site_id',values='wbt').reindex(index=dates,columns=primary.site_id).to_numpy()
 b=d.pivot(index='analysis_date',columns='site_id',values='wbt').reindex(index=dates,columns=grid.dense_site_id).to_numpy()
 assert np.isfinite(a).all() and np.isfinite(b).all()
 err=float(np.max(abs(a-b[:,pidx])))
 assert err<1e-6, (year,err)
 ps=p.groupby('analysis_date').first().reindex(dates);ds=d.groupby('analysis_date').first().reindex(dates)
 assert (ps.time.to_numpy()==ds.time.to_numpy()).all()
 labelerr=float(np.max(abs(ps.regional_mean_wbt-ds.label_mean_wbt)))
 assert labelerr<1e-6
 assert np.max(abs(a.mean(axis=1)-ps.regional_mean_wbt))<1e-5
 days=pd.DataFrame(dict(analysis_date=dates,year=year,month=pd.to_datetime(dates).month,time=ps.time.to_numpy(),regional_mean_wbt=ps.regional_mean_wbt.to_numpy()))
 days['record_id']=days.year*100+days.month
 for month,indices in days.groupby('month').groups.items():
  lo,hi=np.quantile(days.loc[indices,'regional_mean_wbt'],[.25,.75],method='linear')
  x=days.loc[indices,'regional_mean_wbt'];days.loc[indices,'regime']=np.where(x>=hi,'high',np.where(x>=lo,'middle','low'))
  days.loc[indices,'q25']=lo;days.loc[indices,'q75']=hi
 days['role']='development' if year in DEVELOPMENT else 'evaluation'
 manifest.append(days);pfields.append(a);dfields.append(b)
 audits.append(dict(year=year,days=92,primary_sites=121,dense_sites=465,embedded_max_error=err,label_mean_max_error=labelerr))
 print('Audited',year,flush=True)
days=pd.concat(manifest,ignore_index=True);a=np.vstack(pfields);b=np.vstack(dfields)
assert sum(days.role=='evaluation')==3036
prior=pd.read_csv(ROOT/'results/sensitivity_event_manifest.csv')
merged=days.merge(prior[['date','regime','peak_time']],left_on='analysis_date',right_on='date',suffixes=('','_stored'),validate='one_to_one')
assert len(merged)==3220 and (merged.regime==merged.regime_stored).all() and (merged.time==merged.peak_time).all()
days.to_parquet(OUT/'day_manifest.parquet',index=False)
np.savez_compressed(OUT/'audited_fields.npz',primary=a,dense=b,primary_dense_indices=pidx)
# Threshold selection uses development fields only, before any evaluation footprints.
dev=b[days.role=='development'].ravel()
ql,qh=np.quantile(dev,[.01,.995],method='linear');step=.25
lo=math.floor(ql/step)*step;hi=math.ceil(qh/step)*step
thresholds=np.round(np.arange(lo,hi+step/2,step),8).tolist()
conf=dict(frozen_at=datetime.now(timezone.utc).isoformat(),development_years=list(DEVELOPMENT),evaluation_years=[y for y in range(1991,2026) if y not in DEVELOPMENT],thresholds_c=thresholds,threshold_step_c=step,threshold_source_quantiles=[.01,.995],development_quantiles_c=[ql,qh],threshold_endpoints='outward rounding; strict exceedance >',domain=[105,20,125,42],land_source='R spData Natural Earth low-resolution world land union, all countries within original rectangle',land_source_sha256=sha256(boundary),adjacency='rook primary; queen sensitivity',area_sensitivity='cosine latitude on retained grid; identical labels and times',threshold_sensitivity='every second threshold at both offsets; no refitting',regression='weighted ridge least squares for proportion responses; predictions clipped to [0,1]',splines=dict(n_knots=4,degree=2,knots='uniform within training range',extrapolation='linear'),ridge_alphas=[.1,1.,10.,100.,1000.],tuning='nested five-year blocked CV with one-calendar-year buffer, same folds all summaries; choose alpha separately per outcome by equal-year/equal-month MAE averaged across all thresholds and all days; one alpha per outer fold/method/outcome',features=dict(mean='primary regional mean and day-of-summer quadratic B-splines',simple='mean + spatial variance + signed OLS latitude slope + primary Voronoi-area weighted exceedance proportion, each with same spline basis',variogram='simple + log(binned semivariance/V), five pair-distance quantile bins fixed by coordinates',graph='simple + log(Q_h/V), original five physical bandwidths'),zero_rule='if V==0, all shape ratios take continuous constant-field convention 0 on log scale; if V>0 any nonpositive Q/bin energy triggers diagnostic failure, no epsilon',interpolation='piecewise linear in original projected kilometre coordinates; nearest original node outside convex hull; no tuning',uncertainty='paired calendar moving-block bootstrap of annual summaries, length 5, 1999 replicates, seed 20260907; pointwise 95% bands; shared resampling all thresholds/outcomes; no daily iid inference; model comparison intervals conditional on fitted cross-validation predictions',models=['mean','simple','variogram','graph','interpolation'],outcomes=['A','C','A_new344'],no_new_simulations=True)
save_json('config.json',conf);save_json('input_hashes.json',inputs)
save_json('data_audit.json',dict(status='passed',years=audits,evaluation_days=3036,development_days=184,embedded_sites=121,new_sites=344,stored_labels_and_times_match=True,land_area_km2=float(grid.land_area_km2.sum()),zero_land_cells=grid.loc[grid.land_area_km2==0,'dense_site_id'].tolist(),rectangle_land_km2=area(land),coverage_fraction_of_rectangle_land=float(grid.land_area_km2.sum()/area(land)),threshold_range_c=[lo,hi],n_thresholds=len(thresholds),original_hash_files=len(protected),config_sha256=sha256(OUT/'config.json')))
print(json.dumps(json.loads((OUT/'data_audit.json').read_text()) | {'years':'35 passed'},indent=2))
