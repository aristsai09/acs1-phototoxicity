from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    balanced_accuracy_score,
)
from scipy.stats import spearmanr, binomtest

ROOT = Path(__file__).resolve().parents[2]
PRED = ROOT / "data" / "frozen_predictions" / "19_core_700tree_compound_mean_oof_predictions.csv"
STATE = ROOT / "data" / "state_tables" / "005_04_stda_state_summary.csv"
OUT = ROOT / "reproduced_results"
OUT.mkdir(exist_ok=True)

pred = pd.read_csv(PRED)
models = {
    "A": "A_structure_only",
    "C": "C_parent_11",
    "D": "D_microstate_summary_56",
    "E": "E_microstate_conformer_184",
    "Dshuf": "D_shuffled_microstate_56",
    "Eshuf": "E_shuffled_microstate_184",
}
y = pred["y"].astype(int).to_numpy()
n = len(y)

rows = []
for short, col in models.items():
    p = pred[col].to_numpy(float)
    yh = (p >= 0.5).astype(int)
    rows.append(
        dict(
            model=short,
            feature_set=col,
            n=n,
            roc_auc=roc_auc_score(y, p),
            pr_auc=average_precision_score(y, p),
            brier=brier_score_loss(y, p),
            f1=f1_score(y, yh),
            mcc=matthews_corrcoef(y, yh),
            balanced_accuracy=balanced_accuracy_score(y, yh),
        )
    )
pd.DataFrame(rows).to_csv(OUT / "compound_mean_metrics.csv", index=False)

B = 10000
rng = np.random.default_rng(20260923)
pos = np.flatnonzero(y == 1)
neg = np.flatnonzero(y == 0)
wp = rng.multinomial(len(pos), np.full(len(pos), 1 / len(pos)), size=B).astype(np.float32)
wn = rng.multinomial(len(neg), np.full(len(neg), 1 / len(neg)), size=B).astype(np.float32)
W = np.zeros((B, n), dtype=np.float32)
W[:, pos] = wp
W[:, neg] = wn
Ptot = W[:, pos].sum(1)
Ntot = W[:, neg].sum(1)
Tot = Ptot + Ntot

def weighted_metrics(prob):
    prob = np.asarray(prob, float)
    err2 = (prob - y) ** 2
    brier = (W @ err2) / Tot

    order = np.argsort(prob, kind="mergesort")
    yy = y[order]
    ww = W[:, order]
    auc_num = np.zeros(B, float)
    cum_neg = np.zeros(B, float)
    vals = prob[order]
    starts = np.r_[0, np.flatnonzero(np.diff(vals) != 0) + 1]
    ends = np.r_[starts[1:], len(order)]
    for s, e in zip(starts, ends):
        wg = ww[:, s:e]
        yg = yy[s:e]
        posw = wg[:, yg == 1].sum(1) if np.any(yg == 1) else 0.0
        negw = wg[:, yg == 0].sum(1) if np.any(yg == 0) else 0.0
        auc_num += posw * (cum_neg + 0.5 * negw)
        cum_neg += negw
    auc = auc_num / (Ptot * Ntot)

    order = np.argsort(-prob, kind="mergesort")
    yy = y[order]
    ww = W[:, order]
    vals = prob[order]
    starts = np.r_[0, np.flatnonzero(np.diff(vals) != 0) + 1]
    ends = np.r_[starts[1:], len(order)]
    cum_p = np.zeros(B, float)
    cum_n = np.zeros(B, float)
    ap = np.zeros(B, float)
    for s, e in zip(starts, ends):
        wg = ww[:, s:e]
        yg = yy[s:e]
        pw = wg[:, yg == 1].sum(1) if np.any(yg == 1) else 0.0
        nw = wg[:, yg == 0].sum(1) if np.any(yg == 0) else 0.0
        cum_p += pw
        cum_n += nw
        precision = np.divide(cum_p, cum_p + cum_n, out=np.zeros_like(cum_p), where=(cum_p + cum_n) > 0)
        ap += precision * (pw / Ptot)
    return {"roc_auc": auc, "pr_auc": ap, "brier": brier}

boot_by_model = {k: weighted_metrics(pred[v].to_numpy(float)) for k, v in models.items()}
comparisons = [
    ("C-A", "A", "C"),
    ("D-C", "C", "D"),
    ("E-D", "D", "E"),
    ("D-A", "A", "D"),
    ("E-A", "A", "E"),
    ("D-Dshuf", "Dshuf", "D"),
    ("E-Eshuf", "Eshuf", "E"),
]
boot_summ = []
for label, frm, to in comparisons:
    for metric in ["roc_auc", "pr_auc", "brier"]:
        draws = boot_by_model[to][metric] - boot_by_model[frm][metric]
        p0 = pred[models[frm]].to_numpy(float)
        p1 = pred[models[to]].to_numpy(float)
        if metric == "roc_auc":
            obs = roc_auc_score(y, p1) - roc_auc_score(y, p0)
            supportive = np.mean(draws > 0)
        elif metric == "pr_auc":
            obs = average_precision_score(y, p1) - average_precision_score(y, p0)
            supportive = np.mean(draws > 0)
        else:
            obs = brier_score_loss(y, p1) - brier_score_loss(y, p0)
            supportive = np.mean(draws < 0)
        lo, hi = np.quantile(draws, [0.025, 0.975])
        boot_summ.append(
            dict(
                comparison=label,
                metric=metric,
                observed=float(obs),
                ci_low=float(lo),
                ci_high=float(hi),
                supportive_fraction=float(supportive),
                n_boot=B,
            )
        )
pd.DataFrame(boot_summ).to_csv(OUT / "stepwise_bootstrap_10000.csv", index=False)

trans = []
for label, frm, to in [("C-A", "A", "C"), ("D-C", "C", "D"), ("E-D", "D", "E"), ("D-A", "A", "D"), ("E-A", "A", "E")]:
    a = (pred[models[frm]].to_numpy(float) >= 0.5).astype(int)
    b = (pred[models[to]].to_numpy(float) >= 0.5).astype(int)
    ca = a == y
    cb = b == y
    rescued = int(np.sum((~ca) & cb))
    harmed = int(np.sum(ca & (~cb)))
    nd = rescued + harmed
    pval = float(binomtest(rescued, n=nd, p=0.5, alternative="two-sided").pvalue) if nd else np.nan
    trans.append(dict(comparison=label, rescued=rescued, harmed=harmed, net_rescue=rescued - harmed, discordant=nd, mcnemar_exact_p=pval))
pd.DataFrame(trans).to_csv(OUT / "classification_transitions.csv", index=False)

st = pd.read_csv(STATE)
st = st[(st["state_trusted"] == 1) & (st["status"].astype(str).str.lower() == "complete")].copy()
features = [
    "s1_eV",
    "t1_eV",
    "first_bright_nm",
    "max_f",
    "sum_f_uv_290_400",
    "max_f_uv_290_400",
    "sum_f_photo_290_700",
    "max_f_photo_290_700",
    "s1_t1_gap_eV",
    "bright_t1_gap_eV",
]
for c in features:
    st[c] = pd.to_numeric(st[c], errors="coerce")
agg = []
for cid, g in st.groupby("compound_id"):
    rec = {"compound_id": cid, "trusted_complete_states": len(g)}
    for c in features:
        v = g[c].dropna().to_numpy(float)
        rec[c + "_range"] = float(np.ptp(v)) if len(v) >= 2 else np.nan
        rec[c + "_sd"] = float(np.std(v, ddof=1)) if len(v) >= 2 else np.nan
    agg.append(rec)
agg = pd.DataFrame(agg)
ana = pred[["compound_id", "name", "y", "C_parent_11", "D_microstate_summary_56", "E_microstate_conformer_184"]].merge(agg, on="compound_id", how="left")
ana["D_minus_C"] = ana["D_microstate_summary_56"] - ana["C_parent_11"]
ana["abs_D_minus_C"] = ana["D_minus_C"].abs()
ana["abs_E_minus_D"] = (ana["E_microstate_conformer_184"] - ana["D_microstate_summary_56"]).abs()
ana.to_csv(OUT / "microstate_variability_master.csv", index=False)

corr_feats = [
    "trusted_complete_states",
    "s1_eV_range",
    "t1_eV_range",
    "first_bright_nm_range",
    "max_f_range",
    "sum_f_uv_290_400_range",
    "max_f_uv_290_400_range",
    "sum_f_photo_290_700_range",
    "max_f_photo_290_700_range",
    "s1_t1_gap_eV_range",
]
cr = []
for c in corr_feats:
    d = ana.loc[ana["trusted_complete_states"] >= 2, [c, "abs_D_minus_C"]].dropna()
    if len(d) >= 10:
        rho, pval = spearmanr(d[c], d["abs_D_minus_C"])
        cr.append(dict(feature=c, n=len(d), spearman_rho=float(rho), p_value=float(pval)))
cr = pd.DataFrame(cr)
if len(cr):
    m = len(cr)
    order = np.argsort(cr.p_value.values)
    q = np.empty(m)
    prev = 1.0
    for ri in range(m - 1, -1, -1):
        i = order[ri]
        rank = ri + 1
        val = min(prev, cr.p_value.iloc[i] * m / rank)
        q[i] = val
        prev = val
    cr["q_BH"] = q
    cr["significant_fdr_0_05"] = cr.q_BH < 0.05
    cr = cr.sort_values("p_value")
cr.to_csv(OUT / "microstate_variability_correlations_authoritative.csv", index=False)

print(f"Wrote reproduced results to {OUT}")
print(pd.DataFrame(rows).to_string(index=False))
