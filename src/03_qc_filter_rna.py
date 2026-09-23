"""Apply snRNA QC thresholds decided in docs/02 and flag doublets per sample.

Rules (per sample unless stated):
  low counts : log1p(total_counts) < median - 5 MAD
  low genes  : log1p(n_genes)      < median - 5 MAD
  high mt    : pct_counts_mt > 5 % (fixed; nuclei should be near 0)
  doublet    : scrublet, run per sample
No upper bound on counts or genes: the high-count mode is neurons (docs/02).
"""
import numpy as np
import pandas as pd
import scanpy as sc

IN = "data/processed/rna_loaded.h5ad"
OUT = "data/processed/rna_qc.h5ad"
TABLE = "results/03_qc_rna_filter_summary.tsv"
MAD_K = 5
MT_MAX = 5.0


def low_outlier(x, k):
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return x < med - k * mad


def main():
    adata = sc.read_h5ad(IN)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

    obs = adata.obs
    obs["flag_low_counts"] = False
    obs["flag_low_genes"] = False
    for s, idx in obs.groupby("sample", observed=True).groups.items():
        obs.loc[idx, "flag_low_counts"] = low_outlier(obs.loc[idx, "log1p_total_counts"].values, MAD_K)
        obs.loc[idx, "flag_low_genes"] = low_outlier(obs.loc[idx, "log1p_n_genes_by_counts"].values, MAD_K)
    obs["flag_high_mt"] = obs["pct_counts_mt"] > MT_MAX

    # scrublet per sample: doublet rate depends on each sample's loading, so
    # each sample gets its own simulated-doublet model and its own threshold
    score = pd.Series(np.nan, index=obs.index)
    pred = pd.Series(False, index=obs.index)
    for s in obs["sample"].unique():
        sub = adata[obs["sample"] == s].copy()
        sc.pp.scrublet(sub, random_state=0, verbose=False)
        score[sub.obs_names] = sub.obs["doublet_score"].values
        pred[sub.obs_names] = sub.obs["predicted_doublet"].astype(bool).values
    obs["doublet_score"] = score
    obs["flag_doublet"] = pred

    flags = ["flag_low_counts", "flag_low_genes", "flag_high_mt", "flag_doublet"]
    obs["qc_pass"] = ~obs[flags].any(axis=1)

    summary = obs.groupby("sample", observed=True).agg(
        n=("qc_pass", "size"),
        low_counts=("flag_low_counts", "sum"),
        low_genes=("flag_low_genes", "sum"),
        high_mt=("flag_high_mt", "sum"),
        doublet=("flag_doublet", "sum"),
        removed=("qc_pass", lambda x: int((~x).sum())),
        kept=("qc_pass", "sum"),
    )
    summary["pct_removed"] = (100 * summary["removed"] / summary["n"]).round(1)
    summary.to_csv(TABLE, sep="\t")
    print(summary.to_string())
    print(f"\ntotal kept {obs['qc_pass'].sum()} / {len(obs)}")

    # agreement with the authors' own QC
    ct = pd.crosstab(obs["author_retained"], obs["qc_pass"], rownames=["author_retained"], colnames=["our_qc_pass"])
    print("\n" + ct.to_string())

    kept = adata[obs["qc_pass"]].copy()
    kept.layers["counts"] = kept.X.copy()
    kept.write_h5ad(OUT)
    print(f"\nwrote {OUT}: {kept.shape}")


if __name__ == "__main__":
    main()
