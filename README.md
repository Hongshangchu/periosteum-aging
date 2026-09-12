# Periosteum aging: vascular–quiescent program loss and abortive stromal state transition

Analysis code accompanying the manuscript:

> **Vascular–quiescent program loss and abortive stromal state transition define the aged periosteum.**
> Hongshang Chu, Jianguo Qi, Yinxuan Guo, Ziyu Liu, Yong Liu.
> Department of Traumatic Orthopedics, Affiliated Hospital of Shandong Second Medical University, Weifang, China.

This repository contains the complete analysis pipeline for the single-cell re-analysis
of young and aged mouse periosteum, its cross-cohort validation, human genetic
fine-mapping (cis-Mendelian randomization and colocalization), and figure generation.
Pseudorandom seeds are fixed inside the scripts.

---

## Repository contents

| Script | Role | Key outputs (manuscript items) |
|---|---|---|
| `run_periosteum_aging20.R` (v3.8) | Main pipeline: preprocessing (Seurat/Harmony), aging DE, GO/GSEA, senescence signatures, pseudotime trajectory, CellChat communication, external-cohort validation, cis-MR + colocalization | Table 1, Table 2; Figures 1–5, Figures S1–S4; Supplementary Tables S1–S23, S25 |
| `run_periosteum_strengthening20-7.R` (v1.1) | Robustness layer: pseudo-library bootstrap, label permutations, multi-signature concordance, cross-method trajectory validation (Slingshot), ligand-activity scoring | Robustness statistics quoted in the text; Supplementary Tables (S12–S19-related assets) |
| `run_coloc_headline4_eBMD.R` | Colocalization (coloc.abf; coloc.susie attempted) of the four Bonferroni-surviving MR genes (DOCK9, SPP1, TMTC2, APLNR) against heel eBMD (primary) and fractures (sensitivity), GTEx v8 fibroblasts / whole blood cis-eQTL | Supplementary Tables S26, S27; Figure 5D source data |
| `periosteum_figures_fixed_v3.6_integrated17-1.R` | Figure assembly for main and supplementary figures | `figures_fixed/*.pdf/.tiff` (Figures 1–5, S1–S5) |

## Data

All input data are publicly available; no new sequencing data were generated.

- **GEO (scRNA-seq / snRNA-seq):**
  - GSE280914 — discovery mouse periosteum atlas (young 3–4 mo vs aged 20–24 mo; intact vs day-3 fracture; CD45⁺/CD45⁻; 8 libraries)
  - GSE297256 — aging bone/cartilage single-cell atlas (external validation)
  - GSE198666 — young/old fracture callus (external validation)
  - GSE232516 — human periosteum single-cell dataset
  - GSE278165 — human periosteal stem/progenitor cell reference
- **OpenGWAS outcome summary statistics:** heel BMD `ebi-a-GCST90029004`, femoral-neck BMD `ieu-a-980`, osteoporosis `ebi-a-GCST90038656`, fractures `ebi-a-GCST90038703`, `ebi-a-GCST90038705`
- **cis-eQTL:** GTEx v8 (Cells_Cultured_fibroblasts, Whole_Blood, Skeletal_Muscle) via TwoSampleMR / direct download from the GTEx v8 public bucket

## Environment

- R 4.4.x (developed and run under R 4.4.2, Windows 11 x64)
- Key packages: Seurat, Harmony, Monocle3, Slingshot, CellChat, clusterProfiler, msigdbr,
  AUCell, UCell, TwoSampleMR, coloc (5.2.3), ieugwasr, susieR, org.Hs.eg.db,
  data.table (1.18.x), tidyverse (2.0.0)
- Hardware notes: extracting GTEx v8 variant annotations for colocalization requires
  ~8–16 GB free RAM (see `run_coloc_headline4_eBMD.R` header)

## Reproduction workflow

Run from a single project root directory (scripts create `output/`, `figures_fixed/`,
`coloc_data/` as needed):

```bash
# 1. Main pipeline (downloads GEO data on first run; longest step)
Rscript run_periosteum_aging20.R

# 2. Robustness layer (depends on step 1 outputs)
Rscript run_periosteum_strengthening20-7.R

# 3. Colocalization of the four Bonferroni-surviving genes vs eBMD/fractures
#    (reuses GTEx cache if present; otherwise downloads GTEx v8 files)
Rscript run_coloc_headline4_eBMD.R

# 4. Figure assembly
Rscript periosteum_figures_fixed_v3.6_integrated17-1.R
```

## Notes

- The discovery cohort provides one library per group; donor-level statistics were
  therefore not possible, and cell-level p-values are flagged as exploratory in the
  manuscript. The 50 high-confidence aging-associated genes (Table S24) were selected
  as the top-ranking aging DEGs (adjusted p, |log2FC| > 0.25) in an exploratory
  cell-level analysis, because sample-level age × injury interaction modelling was
  not feasible with fewer than two libraries per stratum.
- coloc.susie could not be computed for any gene–tissue pair (in-sample LD reference
  unavailable); coloc.abf results are reported in Tables S26–S27.
- The GTEx v8 eQTL extraction reuses/downloads cached files under `coloc_data/_cache/`;
  these are large and are not part of this repository.

## Citation

If you use this code, please cite the manuscript (above). A citable DOI for this
repository will be added upon publication.

## License

[MIT](LICENSE)
