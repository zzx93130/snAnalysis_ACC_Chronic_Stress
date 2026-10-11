"""Step 7 - balanced cell sets for CellChat (logic of script/cellchat_split3.R).

Cell groups with < 30 cells in any biological sample are excluded; the remaining groups are
downsampled to the same number of cells in every sample: min(smallest sample count, 500).
Two groupings:
  fine   : neuron subtypes + glial/vascular major types   (Fig 3C, differential count/weight)
  coarse : Cdh9+, OtherExNeuron, InNeuron + glial/vascular major types   (Fig 2G)
Writes data/cellchat_<cfg>/{matrix.mtx, genes.txt, meta.csv, balancing.csv}
(matrix = log-normalised expression, genes x cells, CellChat's documented input).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.io as sio
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA, SEED
from config import *  # project paths

MIN_CELLS, MAX_CELLS = 30, 500
adata = sc.read_h5ad(H5AD / "ACC_all_cells.h5ad")
del adata.layers["counts"]
obs = adata.obs
coarse = obs["subtype"].astype(str).where(obs["subtype"].astype(str) == "Cdh9+",
                                          obs["major_celltype"].astype(str).replace({"ExNeuron": "OtherExNeuron"}))
groupings = {"fine": obs["subtype"].astype(str), "coarse": coarse}
rng = np.random.default_rng(SEED)
for cfg, lab in groupings.items():
    tab = pd.crosstab(lab, obs["sample"])
    keep = tab.index[(tab >= MIN_CELLS).all(1)]
    target = tab.loc[keep].min(1).clip(upper=MAX_CELLS)
    cells = []
    for ct in keep:
        for s in tab.columns:
            idx = np.where((lab == ct).values & (obs["sample"] == s).values)[0]
            cells.extend(rng.choice(idx, int(target[ct]), replace=False))
    cells = np.sort(np.array(cells))
    out = DATA / f"cellchat_{cfg}"; out.mkdir(exist_ok=True)
    sub = adata[cells]
    sio.mmwrite(out / "matrix.mtx", sub.X.T.tocsc())
    pd.Series(sub.var_names).to_csv(out / "genes.txt", index=False, header=False)
    meta = pd.DataFrame({"cell_group": lab.iloc[cells].values, "sample": sub.obs["sample"].values,
                         "condition": sub.obs["Experiment"].astype(str).values}, index=sub.obs_names)
    meta.to_csv(out / "meta.csv")
    bal = tab.copy(); bal["kept"] = bal.index.isin(keep); bal["per_sample_target"] = target.reindex(bal.index)
    bal.to_csv(out / "balancing.csv")
    print(cfg, len(keep), "groups kept;", "excluded:", list(tab.index.difference(keep)), "cells", len(cells))
