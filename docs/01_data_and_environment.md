# 01 — Data source and compute environment

Recorded 2026-09-22. Everything below was verified by direct inspection on the day it was written.

## Dataset

GSE174367: Morabito et al. 2021, *Nature Genetics*. Single-nucleus RNA-seq and single-nucleus ATAC-seq of
human prefrontal cortex, late-stage Alzheimer's disease vs age-matched controls. The two modalities were
generated from separate nuclei of the same donors, not from the same nuclei (this is not a 10x Multiome
dataset). Cross-modality alignment is therefore computational, at the level of donor and cell type.

Reasons for choosing this dataset over PBMC-type benchmarks: human brain is a tissue where I can judge whether
a result is biologically plausible or an artifact, and it connects to the in vivo mouse hippocampus work in
my dissertation.

## Files on GEO

All five supplementary files are under `ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/`.

| File | Bytes | Status |
|---|---|---|
| `GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` | 273,975,534 | downloaded to `data/raw/`, byte count verified |
| `GSE174367_snRNA-seq_cell_meta.csv.gz` | 435,170 | `data/` |
| `GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5` | 360,317,403 | downloaded to `data/raw/`, byte count verified |
| `GSE174367_snATAC-seq_cell_meta.csv.gz` | 1,066,930 | `data/` |
| `GSE174367_bulkRNA_processed.rda.gz` | — | not downloaded; bulk RNA-seq from the same donors, reserved for a possible deconvolution phase |

`data/raw` is a symlink to `/scratch/xw347/ad-cortex-multiomic/raw`. No data file is tracked in git.

## Cell metadata

Column order differs between the two tables.

- snRNA: `Barcode, SampleID, Diagnosis, Batch, Cell.Type, cluster, Age, Sex, PMI, Tangle.Stage, Plaque.Stage, RIN`
- snATAC: `Sample.ID, Batch, Age, Sex, PMI, Tangle.Stage, Plaque.Stage, Diagnosis, RIN, cluster, Cell.Type, Barcode`

Both tables carry the authors' `Cell.Type` labels. These are used as one of two annotation tracks for
comparison, not as ground truth.

Barcode suffixes differ: snRNA barcodes all end in `-1` and may not be unique across samples; snATAC barcodes
carry a per-sample suffix (e.g. `-13`). Sample assignment therefore relies on the metadata table, not on the
barcode string. Barcode uniqueness inside each h5 is checked in the next step before any QC.

## Experimental design

Counts below are nuclei per sample from the metadata tables (pre-QC).

| | snRNA-seq | snATAC-seq |
|---|---|---|
| Nuclei | 61,472 | 130,418 |
| Donors | 18 (11 AD, 7 control) | 20 (12 AD, 8 control) |
| Batch 1 | 2 AD, 2 control | 2 AD, 2 control |
| Batch 2 | 4 AD, 2 control | 5 AD, 3 control |
| Batch 3 | 5 AD, 3 control | 5 AD, 3 control |
| Nuclei per sample | 2,544 to 4,835 | 2,605 to 11,971 |

Every batch contains both diagnoses in both modalities, so batch and diagnosis are not confounded and can be
modelled together in differential expression (`~ Batch + Diagnosis`).

The two snATAC-only donors are Sample-101 (control) and Sample-40 (AD), both in Batch 2.

In snRNA, the four Batch 1 samples each have 4,400 to 4,800 nuclei while the other two batches mostly fall
between 2,500 and 3,900. In snATAC the spread is wider (2,605 to 11,971). Per-sample nuclei counts differ
enough that QC thresholds will be set per sample rather than globally.

## Constraints imposed by the released files

1. **Filtered matrices only, no raw matrix.** CellBender needs the raw (unfiltered) matrix to model ambient
   RNA and cannot be used here. Ambient correction will use a method that works on the filtered matrix
   (decontX or SoupX in filtered mode); its estimates are less precise. Re-processing from SRA FASTQ with
   Cell Ranger would recover the raw matrix at a cost of one to two weeks and several hundred GB of scratch.
2. **snATAC released as a peak-by-cell matrix, no fragments file.** Peaks are the authors' definitions and
   cannot be re-called. Fragment-level QC metrics (TSS enrichment, nucleosome signal) cannot be computed.
   Peak-matrix-level QC and MultiVI integration are unaffected.
3. **Each h5 is a single file aggregated across samples.** Per-sample splitting depends on the metadata
   table (see Cell metadata above).

## Decisions

- **Analyse the 18 donors present in both modalities.** The multimodal integration can only use the
  intersection, and keeping the same donor set everywhere means one sample table, one pseudobulk `n`, and
  one description in every document. Sample-101 and Sample-40 are excluded from the snATAC analysis after
  this evaluation, not for any quality reason.
- **Start from the released filtered matrices.** The full workflow is built and documented on these first.
  The FASTQ route is recorded under "Not yet done" and revisited only if ambient contamination turns out to
  matter for the conclusions.
- **Bulk RNA-seq is out of scope for this phase.** It is the natural second phase (snRNA reference
  deconvolution of the same donors' bulk profiles) and is noted so the option is not forgotten.

## Compute environment

Rutgers Amarel (SLURM). Code lives in `~/ad-cortex-multiomic` (home, tracked); data and intermediates in
`/scratch/xw347/ad-cortex-multiomic` (not tracked). The conda environment `sc` is pinned in
`environment.yml`; versions are what the solver produced on 2026-09-22, frozen after the GPU import test
below passed.

GPU partition: A100-PCIE-40GB nodes, 2 to 4 cards each, NVIDIA driver 570.148.08, RHEL 9.6. Three
node-level issues were found while validating the environment and are handled in `workflow/gpu.slurm`:

| Issue | Symptom | Cause | Handling |
|---|---|---|---|
| PyTorch CUDA 13 wheel | `import torch` hangs without error | CUDA 13 runtime needs driver 580+, nodes have 570 | install the `cu128` wheel (`--extra-index-url https://download.pytorch.org/whl/cu128`) |
| System `libstdc++` | `ImportError: /lib64/libstdc++.so.6: version CXXABI_1.3.15 not found` from scipy | when torch is imported first it loads the system C++ runtime, which is older than what conda-forge scipy needs | `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH` before Python starts |
| `gpuk*` nodes | `import torch` exceeds 90 s, stack in `inspect.getmodule` / `realpath` | these nodes read the home filesystem slowly; torch import touches thousands of files | `--exclude=gpuk[001-018]` |

Some GPU nodes do not isolate cards at the cgroup level; SLURM sets `CUDA_VISIBLE_DEVICES` instead. The
template echoes it so the log shows exactly one device was assigned.

Validation (job 61762915, node gpu015): `torch.cuda.is_available()` returned `True` and `import scvi`
succeeded. The full import chain took about 3 minutes on the GPU node, a fixed startup cost per job.

## Not yet done

- Re-process from SRA FASTQ with Cell Ranger to obtain raw matrices for CellBender. Cost: one to two weeks,
  several hundred GB of scratch. Only worth it if ambient contamination affects conclusions.
- Bulk RNA-seq deconvolution using the snRNA reference (`GSE174367_bulkRNA_processed.rda.gz`).
- Moving the conda environment or a container image to scratch if the 3-minute import cost becomes a
  problem for iteration.
