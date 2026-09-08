"""Fixed, nested blocked comparison of 121-site summaries with 465-site coverage."""
import importlib.util,json,time
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay,cKDTree
from sklearn.preprocessing import SplineTransformer,StandardScaler
from common import ROOT,OUT,CODE,config,load,save_json,hierarchy_weights,verify_originals,sha256
spec=importlib.util.spec_from_file_location('footprints',CODE/'01_footprints.py');fp=importlib.util.module_from_spec(spec);spec.loader.exec_module(fp)
METHODS=['mean','simple','variogram','graph']
ENDPOINTS=['A','C','A_new344']

def geometry(primary,cells):
    lat0=np.deg2rad(primary.lat.mean())
    xy=lambda f:np.column_stack([f.lon.to_numpy()*111.32*np.cos(lat0),f.lat.to_numpy()*110.57])
    x=xy(primary);target=xy(cells)
    dist=np.sqrt(((x[:,None,:]-x[None,:,:])**2).sum(axis=2));upper=np.triu_indices(len(x),1)
    return x,target,dist,upper

def prepare(days,cells,p):
    primary=pd.read_csv(ROOT/'data/grid/eastern_china_121_sites.csv').sort_values('site_id')
    x,target,dist,upper=geometry(primary,cells)
    u=np.array(config()['thresholds_c']);v=p.var(axis=1)
    lat=primary.lat.to_numpy()-primary.lat.mean();grad=(p@lat)/(lat@lat)
    ds=dist[upper];h=np.median(ds)*np.array([.125,.25,.5,1,2])
    sq=(p[:,upper[0]]-p[:,upper[1]])**2
    W=np.exp(-ds[:,None]**2/(2*h[None,:]**2));Q=sq@W/(2*W.sum(axis=0))
    bins=np.quantile(ds,np.arange(1,5)/5);b=np.searchsorted(bins,ds,side='right')
    gamma=np.column_stack([sq[:,b==k].mean(axis=1)/2 for k in range(5)])
    logq=np.zeros_like(Q);logg=np.zeros_like(Q);valid=v>0
    assert np.all(Q[valid]>0) and np.all(gamma[valid]>0)
    logq[valid]=np.log(Q[valid]/v[valid,None]);logg[valid]=np.log(gamma[valid]/v[valid,None])
    near=cKDTree(x).query(target)[1];pw=np.bincount(near,weights=cells.weight,minlength=121)
    direct=np.einsum('i,diu->du',pw,p[:,:,None]>u[None,None,:])
    pd.DataFrame(dict(site_id=primary.site_id,representative_weight=pw)).to_csv(OUT/'primary_cell_weights.csv',index=False)
    feature=days[['analysis_date','year','month','regime']].copy();feature['variance']=v;feature['signed_latitude_gradient']=grad
    for k in range(5):feature[f'Q_{k+1}']=Q[:,k];feature[f'gamma_{k+1}']=gamma[:,k];feature[f'graph_shape_{k+1}']=logq[:,k];feature[f'variogram_shape_{k+1}']=logg[:,k]
    feature.to_parquet(OUT/'primary_features.parquet',index=False)
    # Reproduce the original summary without changing any original results.
    temp=days.copy();temp[[f'q{k}' for k in range(5)]]=Q
    grouped=temp[(temp.role=='evaluation')&temp.regime.isin(['high','middle'])].groupby(['year','month','regime'])[[f'q{k}' for k in range(5)]].mean().unstack('regime')
    effects=np.column_stack([grouped[f'q{k}']['high']/grouped[f'q{k}']['middle']-1 for k in range(5)])
    original=pd.read_csv(ROOT/'results/confirmatory_primary_results.csv');expected=original[original.day_definition=='utc'].estimate.iloc[0]
    err=abs(effects.mean()-expected);assert err<1e-10
    # Full-field interpolation; deterministic coordinates, no reference values.
    tri=Delaunay(x);simplex=tri.find_simplex(target);inside=simplex>=0
    transform=tri.transform[simplex[inside]];bary=np.einsum('ijk,ik->ij',transform[:,:2],target[inside]-transform[:,2]);bary=np.column_stack([bary,1-bary.sum(axis=1)])
    M=np.zeros((len(target),len(x)));ids=np.where(inside)[0]
    for j in range(3):M[ids,tri.simplices[simplex[inside],j]]=bary[:,j]
    M[np.where(~inside)[0],near[~inside]]=1
    assert np.allclose(M.sum(axis=1),1) and M.min()>-1e-9
    assert np.max(abs((p@M.T)[:,np.load(OUT/'audited_fields.npz')['primary_dense_indices']]-p))<1e-8
    interp=p@M.T
    A,C=fp.footprints(interp,cells,u)
    new=~cells.is_original_site.astype(bool).to_numpy();w=cells.weight.to_numpy()[new];w/=w.sum()
    Anew=np.einsum('i,diu->du',w,interp[:,new,None]>u[None,None,:])
    np.savez_compressed(OUT/'interpolation_predictions.npz',A=A,C=C,A_new344=Anew)
    save_json('feature_audit.json',dict(primary_reproduction_error=err,zero_variance_days=int(sum(~valid)),bandwidths_km=h.tolist(),variogram_bin_edges_km=[float(ds.min()),*bins.tolist(),float(ds.max())],bin_counts=np.bincount(b).tolist(),interpolation_extrapolated_dense_sites=int(sum(~inside)),sole_feature_source='121 primary sites; dense coordinates/area used as fixed geometry',original_files=verify_originals()))
    season=(pd.to_datetime(days.analysis_date)-pd.to_datetime(days.year.astype(str)+'-06-01')).dt.days.to_numpy()
    baseline=np.column_stack([days.regional_mean_wbt,season,v,grad])
    return baseline,direct,logg,logq,np.stack([A,C,Anew],axis=2)

def spline_columns(values,tr):
    parts=[]
    for j in range(values.shape[1]):
        z=values[:,j:j+1]
        if np.ptp(z[tr])==0:continue
        model=SplineTransformer(n_knots=4,degree=2,knots='uniform',include_bias=False,extrapolation='linear')
        model.fit(z[tr]);parts.append(model.transform(z))
    return np.column_stack(parts) if parts else np.zeros((len(values),1))

def ridge_predictions(X,Y,tr,va,days,alphas):
    # All fitted transformations use only this split's training observations.
    scaler=StandardScaler().fit(X[tr]);train=scaler.transform(X[tr]);test=scaler.transform(X[va])
    weights=hierarchy_weights(days.iloc[tr])*len(tr)
    mx=np.average(train,axis=0,weights=weights);my=np.average(Y[tr],axis=0,weights=weights)
    train=train-mx;test=test-mx
    gram=train.T@(weights[:,None]*train);rhs=train.T@(weights[:,None]*(Y[tr]-my))
    eig,U=np.linalg.eigh(gram);eig=np.maximum(eig,0)
    transformed=U.T@rhs
    coefs=U@ (transformed[None,:,:]/(eig[None,:,None]+np.array(alphas)[:,None,None]))
    pred=np.einsum('ij,ajk->aik',test,coefs)+my
    return np.clip(pred,0,1)

def run():
    conf=config();days,cells,p,d=load();start=time.time()
    baseline,direct,lg,lq,interp=prepare(days,cells,p)
    z=np.load(OUT/'footprints.npz');Y=np.stack([z[n] for n in ENDPOINTS],axis=2)
    years=days.year.to_numpy();evaluation=days.role.to_numpy()=='evaluation';u=np.array(conf['thresholds_c']);alphas=conf['ridge_alphas']
    blocks=[(a,a+4) for a in range(1991,2026,5)]
    predictions={m:np.full_like(Y,np.nan) for m in METHODS};predictions['interpolation']=interp
    manifests=[];tuning=[];choices=[];inconsistency=[]
    for outer,(low,high) in enumerate(blocks):
        test=np.where(evaluation&(years>=low)&(years<=high))[0]
        train=np.where(evaluation&((years<low-1)|(years>high+1)))[0]
        assert len(set(years[train])&set(years[test]))==0
        for year in range(1991,2026):
            role='development' if year in conf['development_years'] else 'test' if low<=year<=high else 'buffer' if low-1<=year<=high+1 else 'train'
            manifests.append(dict(outer_fold=outer+1,inner_fold=0,year=year,role=role))
        losses={m:np.zeros((len(alphas),3)) for m in METHODS};denom=0
        for inner,(il,ih) in enumerate(blocks):
            va=train[(years[train]>=il)&(years[train]<=ih)]
            tr=train[(years[train]<il-1)|(years[train]>ih+1)]
            if not len(va):continue
            assert len(np.unique(years[tr]))>=5
            for year in np.unique(years[train]):
                role='validation' if il<=year<=ih else 'buffer' if il-1<=year<=ih+1 else 'train'
                manifests.append(dict(outer_fold=outer+1,inner_fold=inner+1,year=int(year),role=role))
            base=spline_columns(baseline,tr);mean=spline_columns(baseline[:,:2],tr)
            nyrs=len(np.unique(years[va]));weight=hierarchy_weights(days.iloc[va]);denom+=nyrs
            for k in range(len(u)):
                directbasis=spline_columns(direct[:,k:k+1],tr)
                strong=np.column_stack([base,directbasis])
                designs=[mean,strong,np.column_stack([strong,lg]),np.column_stack([strong,lq])]
                for method,X in zip(METHODS,designs):
                    pred=ridge_predictions(X,Y[:,k,:],tr,va,days,alphas)
                    mae=np.einsum('aij,i->aj',abs(pred-Y[va,k,:]),weight)
                    losses[method]+=mae*nyrs/len(u)
            print(f'Outer {outer+1}/7 inner {inner+1}/7 complete ({time.time()-start:.0f}s)',flush=True)
        selected={}
        for method in METHODS:
            score=losses[method]/denom;selected[method]=np.argmin(score,axis=0)
            for j,end in enumerate(ENDPOINTS):
                choices.append(dict(outer_fold=outer+1,method=method,outcome=end,alpha=alphas[selected[method][j]]))
                for ai,alpha in enumerate(alphas):tuning.append(dict(outer_fold=outer+1,method=method,outcome=end,alpha=alpha,inner_mae=float(score[ai,j])))
        base=spline_columns(baseline,train);mean=spline_columns(baseline[:,:2],train)
        for k in range(len(u)):
            directbasis=spline_columns(direct[:,k:k+1],train);strong=np.column_stack([base,directbasis])
            for method,X in zip(METHODS,[mean,strong,np.column_stack([strong,lg]),np.column_stack([strong,lq])]):
                preds=ridge_predictions(X,Y[:,k,:],train,test,days,alphas)
                predictions[method][test,k,:]=np.column_stack([preds[selected[method][j],:,j] for j in range(3)])
        print(f'OUTER {outer+1} finished',flush=True)
        np.savez_compressed(OUT/'comparison_checkpoint.npz',**predictions)
    pd.DataFrame(manifests).to_csv(OUT/'fold_manifest.csv',index=False)
    pd.DataFrame(tuning).to_csv(OUT/'tuning_losses.csv',index=False);pd.DataFrame(choices).to_csv(OUT/'selected_penalties.csv',index=False)
    tables=[]
    ids=np.where(evaluation)[0]
    for method,P in predictions.items():
        assert np.isfinite(P[evaluation]).all() and P[evaluation].min()>=-1e-12 and P[evaluation].max()<=1+1e-12
        f=days.iloc[ids].loc[days.iloc[ids].index.repeat(len(u))].reset_index(drop=True)
        f['threshold_c']=np.tile(u,len(ids));f['method']=method;f['outer_fold']=((f.year-1991)//5+1)
        for j,end in enumerate(ENDPOINTS):f[end+'_observed']=Y[ids,:,j].ravel();f[end+'_predicted']=P[ids,:,j].ravel()
        tables.append(f)
        inconsistency.append(dict(method=method,fraction_C_above_A=float(np.mean(P[ids,:,1]>P[ids,:,0]+1e-10)),max_C_above_A=float(np.maximum(P[ids,:,1]-P[ids,:,0],0).max()),fraction_A_increasing_with_threshold=float(np.mean(np.diff(P[ids,:,0],axis=1)>1e-10)),fraction_C_increasing_with_threshold=float(np.mean(np.diff(P[ids,:,1],axis=1)>1e-10))))
    pd.concat(tables,ignore_index=True).to_parquet(OUT/'heldout_predictions.parquet',index=False)
    np.savez_compressed(OUT/'heldout_prediction_arrays.npz',**predictions)
    pd.DataFrame(inconsistency).to_csv(OUT/'prediction_consistency.csv',index=False)
    save_json('comparison_audit.json',dict(status='passed',evaluation_days=int(evaluation.sum()),thresholds=len(u),outer_folds=7,methods=5,no_missing_predictions=True,elapsed_seconds=time.time()-start,config_sha256=sha256(OUT/'config.json'),protocol_sha256=sha256(OUT/'revision_protocol.md'),original_files=verify_originals()))
if __name__=='__main__':run()
