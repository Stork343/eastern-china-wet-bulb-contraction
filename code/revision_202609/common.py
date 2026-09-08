from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output_revision_application'
CODE = Path(__file__).resolve().parent
DEVELOPMENT = (2015, 2022)

def sha256(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def save_json(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False))

def config():return json.loads((OUT/'config.json').read_text())

def load():
    z=np.load(OUT/'audited_fields.npz')
    return pd.read_parquet(OUT/'day_manifest.parquet'),pd.read_csv(OUT/'cell_weights.csv'),z['primary'],z['dense']

def hierarchy_weights(days):
    n=days.groupby(['year','month']).analysis_date.transform('count').to_numpy()
    w=1/n
    return w/w.sum()

def verify_originals():
    original=json.loads((OUT/'original_hashes.json').read_text())
    changed=[name for name,digest in original.items() if not (ROOT/name).exists() or sha256(ROOT/name)!=digest]
    assert not changed, changed
    return dict(files=len(original),changed=changed)
