"""Load GSE174367 snRNA and snATAC matrices, attach sample metadata, write h5ad.

Keeps every nucleus in the released filtered h5 files. Sample assignment comes from
the barcode suffix (verified one-to-one with SampleID in the metadata). Nuclei that
the authors dropped in their own QC are kept and flagged (author_retained = False)
so our QC can be compared against theirs. Restricts both modalities to the 18 donors
present in both.
"""
import pandas as pd
import scanpy as sc

RAW = "data/raw"
META = "data"
OUT = "data/processed"

MODALITIES = {
    "rna": ("GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
            "GSE174367_snRNA-seq_cell_meta.csv.gz", "SampleID"),
    "atac": ("GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5",
             "GSE174367_snATAC-seq_cell_meta.csv.gz", "Sample.ID"),
}
META_COLS = ["Diagnosis", "Batch", "Age", "Sex", "PMI", "Tangle.Stage", "Plaque.Stage", "RIN"]


def load(mod):
    h5, meta_csv, sample_col = MODALITIES[mod]
    adata = sc.read_10x_h5(f"{RAW}/{h5}", gex_only=False)
    adata.var_names_make_unique()
    meta = pd.read_csv(f"{META}/{meta_csv}")
    meta = meta.rename(columns={sample_col: "sample"})
    meta["suffix"] = meta["Barcode"].str.split("-").str[1]

    # suffix -> sample must be one-to-one, otherwise assignment is ambiguous
    per_suffix = meta.groupby("suffix")["sample"].nunique()
    assert (per_suffix == 1).all(), f"{mod}: suffix maps to >1 sample"
    suffix_to_sample = meta.drop_duplicates("suffix").set_index("suffix")["sample"]

    adata.obs["suffix"] = adata.obs_names.str.split("-").str[1]
    adata.obs["sample"] = adata.obs["suffix"].map(suffix_to_sample)
    assert adata.obs["sample"].notna().all(), f"{mod}: unmapped suffix"

    # donor-level covariates, one row per sample
    donor = meta.drop_duplicates("sample").set_index("sample")[META_COLS]
    for col in META_COLS:
        adata.obs[col] = adata.obs["sample"].map(donor[col]).values

    adata.obs["author_retained"] = adata.obs_names.isin(meta["Barcode"])
    author_ct = meta.set_index("Barcode")["Cell.Type"]
    adata.obs["author_celltype"] = adata.obs_names.map(author_ct).fillna("not_retained").values
    adata.obs["modality"] = mod
    return adata


def main():
    rna = load("rna")
    atac = load("atac")

    shared = sorted(set(rna.obs["sample"]) & set(atac.obs["sample"]))
    print(f"donors in both modalities: {len(shared)}")
    for name, adata in (("rna", rna), ("atac", atac)):
        before = adata.n_obs
        adata = adata[adata.obs["sample"].isin(shared)].copy()
        print(f"{name}: {before} -> {adata.n_obs} nuclei after restricting to shared donors; "
              f"{adata.obs['author_retained'].sum()} author-retained")
        adata.write_h5ad(f"{OUT}/{name}_loaded.h5ad")


if __name__ == "__main__":
    main()
