"""Coverage and grid-connected area from synchronous observed fields."""
import numpy as np
import pandas as pd
from scipy.ndimage import label,generate_binary_structure
from common import OUT,config,load,save_json,sha256,verify_originals

def footprints(fields,cells,thresholds,queen=False,weight='weight'):
    w=cells[weight].to_numpy();positive=cells.land_area_km2.to_numpy()>0
    w=w*positive;w=w/w.sum()
    rows,cols=cells.row.to_numpy(),cells.col.to_numpy()
    shape=(rows.max()+1,cols.max()+1)
    A=np.empty((len(fields),len(thresholds)));C=np.empty_like(A)
    structure=generate_binary_structure(2,2 if queen else 1)
    lattice=np.zeros(shape,dtype=bool)
    for i,field in enumerate(fields):
        exceed=(field[:,None]>thresholds[None,:])&positive[:,None]
        A[i]=w@exceed
        for k in range(len(thresholds)):
            lattice[:]=False;lattice[rows,cols]=exceed[:,k]
            components,n=label(lattice,structure)
            masses=np.bincount(components[rows,cols],weights=w,minlength=n+1)
            C[i,k]=masses[1:].max() if n else 0.
    assert np.all(C>=-1e-12) and np.all(C<=A+1e-12) and np.all(A<=1+1e-12)
    assert np.max(np.diff(A,axis=1))<1e-12 and np.max(np.diff(C,axis=1))<1e-12
    return A,C

def long_frame(days,thresholds,arrays):
    n=len(thresholds)
    f=days.loc[days.index.repeat(n)].reset_index(drop=True)
    f['threshold_c']=np.tile(thresholds,len(days))
    for name,a in arrays.items():f[name]=a.ravel()
    return f

if __name__=='__main__':
    conf=config();days,cells,p,d=load();u=np.array(conf['thresholds_c'])
    save_json('checkpoint0.json',dict(status='passed',config_sha256=sha256(OUT/'config.json'),protocol_sha256=sha256(OUT/'revision_protocol.md'),data_audit_sha256=sha256(OUT/'data_audit.json'),original_files=verify_originals()))
    # Meaningful topology test: diagonally separated cells join only under queen adjacency.
    toy=pd.DataFrame(dict(row=[0,0,1,1],col=[0,1,0,1],weight=[.25]*4,land_area_km2=[1]*4))
    a,c=footprints(np.array([[2.,0.,0.,2.]]),toy,np.array([-1.,1.,3.]))
    aq,cq=footprints(np.array([[2.,0.,0.,2.]]),toy,np.array([-1.,1.,3.]),queen=True)
    assert np.allclose(a,[[1,.5,0]]) and np.allclose(c,[[1,.25,0]]) and np.allclose(cq,[[1,.5,0]])
    A,C=footprints(d,cells,u);print('Rook done',flush=True)
    _,C8=footprints(d,cells,u,queen=True);print('Queen done',flush=True)
    Ac,Cc=footprints(d,cells,u,weight='coslat_weight');print('Cosine weighting done',flush=True)
    new=~cells.is_original_site.astype(bool).to_numpy();w=cells.weight.to_numpy()[new];w/=w.sum()
    Anew=np.einsum('i,diu->du',w,d[:,new,None]>u[None,None,:])
    total=cells.land_area_km2.sum()
    arrays=dict(A=A,C=C,S_km2=A*total,largest_km2=C*total,C8=C8,A_coslat=Ac,C_coslat=Cc,A_new344=Anew)
    long_frame(days,u,arrays).to_parquet(OUT/'daily_footprints.parquet',index=False)
    np.savez_compressed(OUT/'footprints.npz',**arrays)
    save_json('footprint_audit.json',dict(status='passed',days=len(days),evaluation_days=int(sum(days.role=='evaluation')),thresholds=len(u),topology_test=True,min_A=float(A.min()),max_A=float(A.max()),max_C_minus_A=float(np.max(C-A)),max_rook_minus_queen=float(np.max(C-C8)),curves_monotone=True,original_files=verify_originals()))
