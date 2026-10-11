"""Step 9 - PROGENy pathway activity per biological sample x cell group
(logic of script/sample_liana_consistent_lr.py: mean log-expression per sample x cell type with
>= 20 cells, ULM with PROGENy mouse top-500 (tmin 5), Welch t-test vs CON per feature,
BH-FDR, feature kept only if present in all samples of CON/CES/CPS).

Outputs: data/progeny_activity_long.csv, data/progeny_differential.csv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, decoupler as dc
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA
from config import *  # project paths

MIN_CELLS = 20
expr = pd.read_csv(DATA / "sample_mean_logexpr_subtype.csv.gz", index_col=0)
meta = pd.read_csv(DATA / "pseudobulk_meta_subtype.csv").set_index("id")
meta = meta[meta.n_cells >= MIN_CELLS]
expr = expr.loc[meta.index]
net = pd.read_csv(DATA / "progeny_mouse_top500.csv.gz")
es, pv = dc.mt.ulm(data=expr, net=net, tmin=5, verbose=False)
act = es.stack().rename("activity").reset_index()
act.columns = ["id", "pathway", "activity"]
act = act.merge(meta[["group", "sample", "condition", "n_cells"]], left_on="id", right_index=True)
# keep group x pathway only if observed in every sample of each condition
nsamp = act.groupby(["group", "condition"])["sample"].nunique().unstack()
ok = nsamp.index[(nsamp[["CON", "CES", "CPS"]] == 3).all(1)]
act = act[act.group.isin(ok)]
act["feature"] = act.group + " | " + act.pathway
act.to_csv(DATA / "progeny_activity_long.csv", index=False)

rows = []
for cond in ["CES", "CPS"]:
    for (g, p), d in act.groupby(["group", "pathway"]):
        a = d.loc[d.condition == "CON", "activity"].values
        b = d.loc[d.condition == cond, "activity"].values
        t = stats.ttest_ind(b, a, equal_var=False)
        rows.append(dict(group=g, pathway=p, feature=f"{g} | {p}", condition=cond,
                         beta=b.mean() - a.mean(), pvalue=t.pvalue))
res = pd.DataFrame(rows)
res["fdr"] = res.groupby("condition")["pvalue"].transform(lambda x: multipletests(x.fillna(1), method="fdr_bh")[1])
w = res.pivot(index=["feature", "group", "pathway"], columns="condition", values=["beta", "pvalue", "fdr"])
w.columns = [f"{a}_{b}" for a, b in w.columns]
w = w.reset_index()
same = np.sign(w.beta_CES) == np.sign(w.beta_CPS)
w["direction"] = np.where(~same, "discordant", np.where(w.beta_CES > 0, "up_in_both", "down_in_both"))
w["effect_distance"] = np.sqrt(w.beta_CES ** 2 + w.beta_CPS ** 2)
w["shared_abs_effect"] = np.minimum(w.beta_CES.abs(), w.beta_CPS.abs())
w = w.sort_values(["effect_distance", "shared_abs_effect"], ascending=False)
w.to_csv(DATA / "progeny_differential.csv", index=False)
print(w.direction.value_counts().to_dict(), "groups", act.group.nunique())
