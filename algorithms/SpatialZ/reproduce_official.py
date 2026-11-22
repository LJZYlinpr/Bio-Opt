import argparse
import hashlib
import inspect
import json
from pathlib import Path
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root / "optimized"))
from cpu import measure, synthesize_cpp
from full_pair import assign_labels
from geometry import load_module
NATIVE = root / "original/SpatialZ.py"

p = argparse.ArgumentParser()
p.add_argument('--starmap', type=Path)
p.add_argument('--out', type=Path, required=True)
opts = p.parse_args()
out = opts.out
out.mkdir(parents=True, exist_ok=True)
if opts.starmap:
    data = ad.read_h5ad(opts.starmap)
    data = data[~data.obs.z.isin(list(range(6,14))+list(range(91,95)))].copy()
    keys = np.array_split(np.unique(data.obs.z), 7)
    frames = [data[data.obs.z.isin(k)].copy() for k in keys]
    assert [a.n_obs for a in frames] == [4073,4187,4169,4102,4110,4162,4175]
    config = {'lr':1e5,'num_projections':200,'Beta':100}
    cases = [(i, i+2, .5, frames[i+1].n_obs) for i in (0,2,4)]
    dataset = 'official real 3D evaluation STARmap'
else:
    directory = root / "small_test/merfish"
    paths = sorted(directory.glob('merfish_*_paste.h5ad'))
    data = ad.concat([ad.read_h5ad(f) for f in paths], join='outer', merge='same')
    keys = ['0.04','0.09','0.14','0.19','0.24']
    frames = [data[data.obs.slice_id.astype(str)==k].copy() for k in keys]
    assert len(paths) == 5 and all(f.n_obs > 0 and f.n_vars == 155 for f in frames)
    config = {'lr':1e-5,'num_projections':80,'Beta':1}
    cases = [(0,1,.75,None),(1,2,.5,None),(3,4,.25,None)]
    dataset = 'official Spatial dense data construction MERFISH'
for f in frames:
    f.obsm['spatial'] = np.asarray(f.obsm['spatial'])[:, :2]
    f.uns.clear()
    for key in list(f.obsp):
        del f.obsp[key]

stock = load_module('stock_spatialz', NATIVE)
fast = load_module('geometry_only_spatialz', NATIVE)
source = NATIVE.read_text()
marker = '    Coor_final = Coor_torch.detach().cpu().numpy()'
assert source.count(marker) == 1
exec(compile(source.replace(marker, marker+'\n    return Coor_final',1), str(NATIVE),'exec'), fast.__dict__)

if config['Beta'] != 100:
    # Reuse the production sampler with the dense tutorial's explicit Beta.
    src = inspect.getsource(synthesize_cpp)
    assert src.count('100 * similarities') == 1
    scope = dict(synthesize_cpp.__globals__)
    scope['demo_beta'] = config['Beta']
    exec(src.replace('100 * similarities', 'demo_beta * similarities', 1), scope)
    synthesize_cpp = scope['synthesize_cpp']

stock_measure = stock.perform_microenv_measure
snap = []
def capture(*args, **kwargs):
    result = stock_measure(*args, **kwargs)
    snap[:] = [np.array(f.obsm['mender'], copy=True) for f in result]
    return result
stock.perform_microenv_measure = capture

torch.set_num_threads(1)
assert torch.ones(2, device='cuda').sum().item() == 2
results = {'dataset': dataset, 'config': {**config,'iterations':1000,'k_sam':50,
           'syn_mode':'default','seed':42,'stock_workers':200,'optimized_workers':32},
           'source_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),
           'gpu':torch.cuda.get_device_name(0), 'cases':[]}
name = 'starmap' if opts.starmap else 'merfish'
for j, (l,r,alpha,n_cell) in enumerate(cases):
    print('CASE START', l,r,alpha,n_cell, flush=True)
    a,b = frames[l],frames[r]
    kwargs=dict(alpha=alpha,n_cell=n_cell,nb_iter_max=1000,seed=42,n_mag=1.,
                syn_mode='default',k_sam=50,device='cuda:0',cell_type_key='cell_class',**config)
    np.random.seed(42)
    t = time.monotonic()
    native = stock.Generate_spatialz(a.copy(),b.copy(),**kwargs)
    stock_time=time.monotonic()-t
    if sp.issparse(native.X):
        native.X=native.X.toarray()
    torch.cuda.empty_cache()
    t=time.monotonic()
    xy = fast.Generate_spatialz(a.copy(),b.copy(),**kwargs)
    labels=assign_labels(a,b,xy)
    query=ad.AnnData(np.zeros(native.shape,dtype=np.float32),
                     obs=pd.DataFrame({'cell_class':pd.Categorical(labels)}),var=a.var.copy())
    query.obsm['spatial']=xy
    aa,bb=a.copy(),b.copy()
    for f in (aa,bb):
        if sp.issparse(f.X):
            f.X=f.X.toarray()
    altered=[aa,bb,query]
    desc=measure(altered,32)
    for f,m in zip(altered,desc,strict=True):
        f.obsm['mender']=m
    arr, _ = synthesize_cpp(altered,seed=42,threads=4)
    fast_time=time.monotonic()-t
    vals=np.asarray(native.X)
    report={'left':l+1,'right':r+1,'alpha':alpha,'cells':len(xy),'genes':len(a.var),
            'xy_exact':bool(np.array_equal(xy,native.obsm['spatial'])),
            'labels_exact':bool(np.array_equal(labels,native.obs.cell_class.astype(str))),
            'mender_exact':[bool(np.array_equal(x,y)) for x,y in zip(snap,desc,strict=True)],
            'expression_exact':bool(np.array_equal(arr,vals)),
            'expression_different_entries':int(np.count_nonzero(arr!=vals)),
            'expression_max_abs_error':float(np.max(np.abs(arr-vals))),
            'stock_wall_seconds':stock_time,'optimized_wall_seconds':fast_time}
    np.savez_compressed(out/f'{name}_{j+1}_comparison.npz', original=vals, optimized=arr,
                        original_xy=native.obsm['spatial'],optimized_xy=xy,
                        original_labels=native.obs.cell_class.to_numpy(dtype=str),optimized_labels=labels)
    results['cases'].append(report)
    dest=out/f'{name}_RESULT.json'
    temp=dest.with_suffix('.tmp')
    temp.write_text(json.dumps(results,indent=2)+'\n')
    temp.replace(dest)
    print(json.dumps(report),flush=True)
