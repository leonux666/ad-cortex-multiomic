"""Integrate across samples with Harmony (baseline) and scVI (main), then score both.

Both methods use sample as the batch variable, not Batch: sample is the finest
unit of technical variation and Batch is nested inside it.
Metrics use the authors' cell-type labels as the biological reference, restricted
to author-retained nuclei so that "not_retained" is not treated as a cell type.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import harmonypy as hm
import numpy as np
import scanpy as sc
import scvi
from scib_metrics.benchmark import Benchmarker

IN = "data/processed/rna_norm.h5ad"
OUT = "data/processed/rna_integrated.h5ad"
MODEL_DIR = "data/processed/scvi_model"
METRICS = "results/05_integration_metrics.tsv"
FIG = "docs/figures/05_umap_integrated.png"
N_LATENT = 30
SEED = 0


def main():
    scvi.settings.seed = SEED
    adata = sc.read_h5ad(IN)
    adata.obs["Batch"] = adata.obs["Batch"].astype(str).astype("category")
    print("input", adata.shape)

    # Harmony: iterative correction of the PCA embedding
    # harmonypy called directly: version 2.0 returns cells x PCs, and scanpy's
    # wrapper still transposes it, so the orientation is checked here instead
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample")
    Z = ho.Z_corr
    Z = Z.cpu().numpy() if hasattr(Z, "cpu") else np.asarray(Z)
    adata.obsm["X_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
    print("harmony done")

    # scVI: trained on raw counts of HVGs, sample as batch covariate
    hvg = adata[:, adata.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(hvg, layer="counts", batch_key="sample")
    model = scvi.model.SCVI(hvg, n_latent=N_LATENT, n_layers=2, gene_likelihood="nb")
    model.train(early_stopping=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    model.save(MODEL_DIR, overwrite=True)
    print("scvi done, epochs trained:", len(model.history["elbo_train"]))

    # UMAPs on each embedding
    for key in ["X_harmony", "X_scVI"]:
        name = key.replace("X_", "")
        sc.pp.neighbors(adata, use_rep=key, key_added=name)
        sc.tl.umap(adata, neighbors_key=name)
        adata.obsm[f"X_umap_{name}"] = adata.obsm["X_umap"].copy()

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    for row, name in enumerate(["harmony", "scVI"]):
        adata.obsm["X_umap"] = adata.obsm[f"X_umap_{name}"]
        for ax, color in zip(axes[row], ["sample", "Batch", "author_celltype"]):
            sc.pl.umap(adata, color=color, ax=ax, show=False, size=3,
                       legend_loc="on data" if color == "author_celltype" else "right margin",
                       legend_fontsize=7, title=f"{name}: {color}")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure", FIG)

    # scib-metrics on author-retained nuclei only
    ret = adata[adata.obs["author_retained"]].copy()
    bm = Benchmarker(ret, batch_key="sample", label_key="author_celltype",
                     embedding_obsm_keys=["X_pca", "X_harmony", "X_scVI"], n_jobs=8)
    bm.benchmark()
    res = bm.get_results(min_max_scale=False)
    res.to_csv(METRICS, sep="\t")
    print(res.round(3).to_string())

    adata.write_h5ad(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
