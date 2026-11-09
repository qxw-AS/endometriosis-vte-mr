# Endometriosis and venous thromboembolism: shared genetic architecture and bidirectional Mendelian randomization

Analysis code for the manuscript *"Mendelian randomization identifies antagonistic pleiotropy at F5 with no detectable genome-wide genetic sharing between endometriosis and venous thromboembolism"* (v1.1, 2026-11-09).

## Overview

The pipeline uses FinnGen release 13 (GRCh38) endometriosis (`N14_ENDOMETRIOSIS`) and venous thromboembolism (`I9_VTE`) summary statistics, with partial external replication in the Rahmioglu et al. (2023) endometriosis meta-analysis (GWAS Catalog `GCST90205183`). It covers:

- LDSC SNP heritability, genetic correlation and the S5 sensitivity analysis excluding F5/F2/ABO windows;
- 58 candidate regions, single-causal-variant colocalization (`coloc.abf`), p12 and sdY prior sensitivity, and PLACO;
- MAGMA gene-set enrichment and a 10,000-permutation gene-level test (seed 20261108);
- bidirectional two-sample MR (IVW, MR-Egger, weighted median, weighted mode, MR-PRESSO) with instrument diagnostics;
- forward and reverse multivariable MR with BMI and ever smoking;
- the F5 rs6025 cross-phenotype lookup, the nine-gene sentinel analysis, external replication and PE/DVT severity stratification.

Primary anchors: rg 0.0443 (95% CI −0.077 to 0.166), F5 PP.H4 0.994, forward MR OR 1.017 (0.975–1.062), reverse MR OR 0.889 (0.840–0.940). Exact values and their manuscript locations are listed in [docs/CODE_MAP.md](docs/CODE_MAP.md).

## Repository layout

```
.
├── README.md
├── config.py                  central path configuration (see below)
├── environment.yml            (not included in this version)
├── 01_sumstats_qc/            download, LDSC munging, HapMap3 SNP tables
├── 02_ldsc_main/              LDSC h2/rg driver — not archived (see Scope)
├── 03_ldsc_sensitivity/       S5 exclusion of F5/F2/ABO +/- 500 kb / 1 Mb
├── 04_regions_coloc/          candidate regions, region extraction, coloc.abf
├── 05_prior_sensitivity/      p12 prior sensitivity at fixed sdY = 0.15
├── 06_placo/                  PLACO, PLACO clumping, integrated summary
├── 07_magma_enrichment/       gene-set enrichment and permutation test
├── 08_instruments_diagnostics/ instrument extraction, harmonisation, diagnostics
├── 09_bidirectional_mr/       bidirectional MR and coagulation-excluded sensitivity
├── 10_mvmr/                   MVMR instrument union, effects, forward and reverse MVMR
├── 11_f5_lookup/              F5 rs6025 across eight FinnGen phenotypes
├── 12_fig3b_validation_severity/ sentinel ORs, external replication, PE/DVT severity
├── 13_endpoints_figures/      endpoint table and all manuscript figures
└── docs/CODE_MAP.md           script -> manuscript section -> expected value
```

## Input data

Place these files in `data/` (or set `ENDO_VTE_DATA`). None are redistributed here.

| File | Source | Content |
|---|---|---|
| `finngen_R13_N14_ENDOMETRIOSIS.gz` | FinnGen R13, https://r13.finngen.fi | 20,534 cases / 122,622 controls |
| `finngen_R13_I9_VTE.gz` | FinnGen R13 | 27,928 cases / 472,258 controls |
| `GCST90205183.h.tsv.gz` (+ `.tbi`) | GWAS Catalog | Rahmioglu 2023 meta-analysis, 60,674 cases / 701,926 controls |
| `1000G_Phase3_ldscores.tgz`, `1000G_Phase3_weights_hm3_no_MHC.tgz`, `1000G_Phase3_plinkfiles.tgz` | Zenodo 10515792 | LDSC and PLINK reference panels |
| `hm3_snps.txt`, `hm3_no_MHC.list.txt` | HapMap3 | SNP lists for munging and MAGMA input |
| `magma_gene_loc.txt`, `ensembl_genes_grch38.tsv`, `GO_BP_2023.gmt` | Ensembl / GO | gene locations and gene sets |
| `finngen_r13_rs6025_api.json` | FinnGen R13 pheweb API (cached) | optional cache for the Fig 3A lookup |

BMI (`BMI_IRN`), ever smoking (`EVER_SMOKER`), pulmonary embolism (`I9_PULMEMB`) and DVT (`I9_PHLETHROMBDVTLOW`) effects are fetched from the FinnGen R13 pheweb API by the MVMR and severity scripts.

## Software

Python 3.11 with numpy, pandas, scipy, matplotlib, intervaltree and requests (`pysam`/tabix optional for remote access in `90_g8_external_validation.py`). External binaries: LDSC 1.0.1, MAGMA 1.10, PLINK 1.9 and liftOver. R 4.3 with TwoSampleMR 0.5.11 and MR-PRESSO 1.0 was used for cross-checks of the MR estimates. Random seed: 20261108.

## Configuration

All paths come from `config.py`:

```python
from config import DATA_DIR, OUT_DIR, FIGURES_DIR
```

Defaults are derived from the location of `config.py`, so the repository can be cloned anywhere. Environment variables override them:

| Variable | Default | Meaning |
|---|---|---|
| `ENDO_VTE_DATA` | `<repo>/data` | summary statistics and reference files |
| `ENDO_VTE_OUT` | `<repo>/out` | analysis outputs |
| `ENDO_VTE_FIGURES` | `<repo>/figures` | Fig 2 and Fig 3 renders |

The four shell scripts in `01_sumstats_qc/` and `07_magma_enrichment/` derive the same root from `BASH_SOURCE` and read `ENDO_VTE_DATA`, so no absolute path needs editing.

## Running

Run the numbered directories in order; later modules read `out/` files written by earlier ones.

```bash
bash 01_sumstats_qc/download_data.sh          # public downloads
bash 01_sumstats_qc/make_ldsc_sumstats.sh     # LDSC-format sumstats
# 02: LDSC h2/rg are run with the ldsc binary; see Scope below
python 03_ldsc_sensitivity/87_s5_ldsc_coag_excl.py
python 04_regions_coloc/80_m3_regions.py
python 04_regions_coloc/81_m3_extract.py
python 04_regions_coloc/82_m3_coloc.py
python 05_prior_sensitivity/86_m3_p12_sens.py
python 06_placo/83_m3_placo.py
python 06_placo/84_m3_placo_clump.py
python 06_placo/85_m3_summary.py
bash 07_magma_enrichment/make_snp_table.sh
bash 07_magma_enrichment/make_magma_pval.sh
python 07_magma_enrichment/12_geneset_perm.py
python 08_instruments_diagnostics/70_mr_extract.py
python 08_instruments_diagnostics/71_mr_harmonize.py
python 08_instruments_diagnostics/73_mr_diagnostics.py
python 09_bidirectional_mr/72_mr_analysis.py
python 09_bidirectional_mr/74_mr_coag_sensitivity.py
python 10_mvmr/96_m5_build_union.py
python 10_mvmr/97_m5_fetch.py
python 10_mvmr/98_m5_mvmr.py
python 10_mvmr/99_m5_reverse_mvmr.py
python 11_f5_lookup/60_fetch_fig3_data.py
python 12_fig3b_validation_severity/62_coag_gene_sentinel.py
python 12_fig3b_validation_severity/90_g8_external_validation.py
python 12_fig3b_validation_severity/99_severity_pe_dvt.py
python 12_fig3b_validation_severity/100_severity_coag_sens.py
python 12_fig3b_validation_severity/114_severity_wmedian.py
python 13_endpoints_figures/113_make_tables.py
python 13_endpoints_figures/20_make_fig2.py
python 13_endpoints_figures/21_make_fig3.py
python 13_endpoints_figures/111_figS1_S2.py
python 13_endpoints_figures/112_figS4_locuszoom.py
python 13_endpoints_figures/115_fig1_workflow.py
```

## Outputs

| Manuscript item | Producer |
|---|---|
| Table I (data sources, LDSC estimates) | `113_make_tables.py`; LDSC values from `out/ldsc/` |
| Table II (instruments) | `113_make_tables.py` from `71_mr_harmonize.py` |
| Table III (MR and MVMR) | `113_make_tables.py` from `72`, `98` |
| Table IV (colocalization and pleiotropy) | `113_make_tables.py` from `82`, `83`, `85` |
| Table V (external replication) | `113_make_tables.py` from `90` |
| Table VI (severity) | `113_make_tables.py` from `99_severity_pe_dvt`, `114` |
| Tables S1A, S1B, S2-S5 | `60`, `62`, `71`, endpoint table, `99_m5_reverse_mvmr` |
| Fig 1 | `115_fig1_workflow.py` (methods-only variant: `116_fig1_methods.py`) |
| Fig 2 | `20_make_fig2.py` |
| Fig 3 | `21_make_fig3.py` |
| Fig S1, S2 | `111_figS1_S2.py` |
| Fig S4 | `112_figS4_locuszoom.py` |

## Scope, known gaps and caveats

- **Module 02 (LDSC h2/rg) has no archived driver script.** LDSC 1.0.1 was invoked directly; the logs (`out/ldsc/h2_endo.log`, `h2_vte.log`, `rg_endo_vte.log`) record the estimates (h2 0.0699 and 0.0185; intercepts 1.139 and 1.121; rg 0.0443).
- **No archived liftOver/QC script.** Build conversion and quality control were run interactively; the munged LDSC inputs are the starting point of the archived pipeline.
- **Not archived:** the three-prevalence-tier analysis, MHC-excluded LDSC, split-half analysis (seed 20261108), the sdY prior sensitivity script, the Risteys/PheWeb endpoint-definition fetch, and the Fig S3 PLACO panel. Their results appear in the manuscript (sdY 0.10-0.30: PP.H4 0.993-0.994; endpoint definitions in Table S3) but cannot be regenerated from this repository alone.
- **Shared controls.** Both primary traits come from the same FinnGen release with a shared control set. The external endometriosis meta-analysis itself contains a FinnGen freeze (6,502 cases), so no SNP-level dataset is FinnGen-free on both sides. Colocalization is therefore reported as within-FinnGen and the external analysis as partial, not independent, replication.

## Data availability

FinnGen release 13: https://r13.finngen.fi. Endometriosis meta-analysis: GWAS Catalog `GCST90205183`. LDSC reference panels: Zenodo record 10515792. Derived tables and figures are regenerated by the scripts above; they are not committed.

## Citation and licence

Citation metadata (`CITATION.cff`) is not included in this version; cite the manuscript above. No licence has been assigned to this repository yet (`LICENSE` is not included).
