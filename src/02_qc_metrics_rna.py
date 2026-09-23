"""Compute per-nucleus QC metrics for snRNA and plot them per sample. No filtering here.

Thresholds are decided after looking at these plots and recorded in docs/02; the
filtering step is a separate script so the decision is visible in the repo history.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

IN = "data/processed/rna_loaded.h5ad"
TABLE = "results/02_qc_rna_per_sample.tsv"
FIG = "docs/figures/02_qc_rna_metrics.png"


def mad_flags(x, k):
    """Return boolean mask of values farther than k MADs from the median."""
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return np.abs(x - med) > k * mad


def main():
    adata = sc.read_h5ad(IN)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True)

    obs = adata.obs
    # per-sample summary: medians and how many nuclei fall outside 5 MADs on log counts / log genes
    rows = []
    for s, g in obs.groupby("sample", observed=True):
        rows.append({
            "sample": s, "batch": g["Batch"].iloc[0], "diagnosis": g["Diagnosis"].iloc[0],
            "n_nuclei": len(g), "author_retained": int(g["author_retained"].sum()),
            "median_counts": int(g["total_counts"].median()),
            "median_genes": int(g["n_genes_by_counts"].median()),
            "median_pct_mt": round(float(g["pct_counts_mt"].median()), 3),
            "p95_pct_mt": round(float(g["pct_counts_mt"].quantile(0.95)), 3),
            "out_counts_5mad": int(mad_flags(g["log1p_total_counts"].values, 5).sum()),
            "out_genes_5mad": int(mad_flags(g["log1p_n_genes_by_counts"].values, 5).sum()),
        })
    table = pd.DataFrame(rows).sort_values(["batch", "diagnosis", "sample"])
    table.to_csv(TABLE, sep="\t", index=False)
    print(table.to_string(index=False))

    # figure: three metrics per sample, samples ordered by batch, coloured by diagnosis
    order = table["sample"].tolist()
    colors = {"AD": "#c0504d", "Control": "#4f81bd"}
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    for ax, metric, label in zip(
        axes,
        ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt"],
        ["log1p(UMI counts)", "log1p(genes detected)", "% mitochondrial"],
    ):
        data = [obs.loc[obs["sample"] == s, metric].values for s in order]
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for body, s in zip(parts["bodies"], order):
            body.set_facecolor(colors[table.set_index("sample").loc[s, "diagnosis"]])
            body.set_alpha(0.7)
        ax.set_ylabel(label)
    axes[-1].set_xticks(range(1, len(order) + 1))
    axes[-1].set_xticklabels([f"{s}\nB{table.set_index('sample').loc[s, 'batch']}" for s in order],
                             rotation=90, fontsize=8)
    axes[0].set_title("snRNA QC metrics per sample (red AD, blue Control; grouped by batch)")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print(f"figure: {FIG}")


if __name__ == "__main__":
    main()
