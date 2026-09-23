"""Normalize, select highly variable genes, build the unintegrated embedding.

Raw counts stay in layers["counts"]; X becomes log1p(library-size-normalized).
HVGs use seurat_v3 on raw counts with batch_key="sample", so a gene has to be
variable within samples to be chosen; genes that differ only between samples are
not favoured. The unintegrated PCA/UMAP is the baseline that 05_integrate is judged against.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

IN = "data/processed/rna_qc.h5ad"
OUT = "data/processed/rna_norm.h5ad"
FIG = "docs/figures/04_umap_unintegrated.png"
N_HVG = 3000
MIN_CELLS = 10
N_PCS = 50


def main():
    adata = sc.read_h5ad(IN)
    print("input", adata.shape)

    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    adata.layers["counts"] = adata.X.copy()
    print("after gene filter", adata.shape)

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_HVG,
                                layer="counts", batch_key="sample")
    print("HVGs", int(adata.var["highly_variable"].sum()))

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata  # full log-normalized matrix kept for marker plots

    sc.pp.pca(adata, n_comps=N_PCS, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_pcs=N_PCS, key_added="unintegrated")
    sc.tl.umap(adata, neighbors_key="unintegrated")
    adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for ax, key in zip(axes, ["sample", "Batch", "author_celltype"]):
        sc.pl.umap(adata, color=key, ax=ax, show=False, size=3, legend_loc="on data" if key == "author_celltype" else "right margin",
                   legend_fontsize=7, title=f"unintegrated UMAP: {key}")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure", FIG)

    adata.write_h5ad(OUT)
    print("wrote", OUT, adata.shape)


if __name__ == "__main__":
    main()
