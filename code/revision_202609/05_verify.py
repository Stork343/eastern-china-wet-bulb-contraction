"""Targeted numerical and separation checks for the application extension."""
import importlib.util,json
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from common import CODE,OUT,load,save_json,hierarchy_weights,verify_originals,sha256
spec=importlib.util.spec_from_file_location('compare',CODE/'03_compare.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
days,cells,p,d=load();Y=np.load(OUT/'footprints.npz')['A'][:,60:63]
X=np.column_stack([p.mean(axis=1),p.var(axis=1),np.arange(len(p))%92])
tr=np.arange(92,276);va=np.arange(30)
ours=m.ridge_predictions(X,Y,tr,va,days,[.1,10.])
s=StandardScaler().fit(X[tr]);weights=hierarchy_weights(days.iloc[tr])*len(tr)
errors=[]
for k,alpha in enumerate([.1,10.]):
 model=Ridge(alpha=alpha).fit(s.transform(X[tr]),Y[tr],sample_weight=weights)
 ref=np.clip(model.predict(s.transform(X[va])),0,1);errors.append(float(abs(ref-ours[k]).max()))
assert max(errors)<1e-9,errors
altered=Y.copy();altered[va]+=10
assert np.allclose(ours,m.ridge_predictions(X,altered,tr,va,days,[.1,10.]),atol=1e-12)
before=json.loads((OUT/'manuscript_before.json').read_text())
from common import ROOT
assert all(sha256(ROOT/name)==digest for name,digest in before.items()),'Manuscript changed before validation finished'
summary=pd.read_csv(OUT/'coverage_summary.csv')
allgrid=summary[summary.outcome.isin(['A','C'])].groupby('outcome').difference.agg(['min','max']);print(allgrid)
save_json('verification.json',dict(status='passed',ridge_matches_sklearn_max_error=max(errors),heldout_response_perturbation_has_no_effect=True,manuscript_unchanged_through_checkpoint3=True,original_files=verify_originals()))
