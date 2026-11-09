#!/usr/bin/env python3
"""Central path configuration for the endometriosis-VTE analysis repository.

No absolute path is hard-coded in the analysis scripts.  Python modules import
the directory constants below; the four shell scripts derive the same root from
their own location and honour the same environment variables.

Default layout (the repository can be cloned anywhere):

    <repo>/data                 input summary statistics and reference files
    <repo>/out                  analysis outputs
    <repo>/out/ldsc             LDSC logs and S5 sensitivity files
    <repo>/out/magma            MAGMA gene and gene-set output
    <repo>/out/m3               candidate regions, coloc, PLACO
    <repo>/out/mr               harmonised instruments, MR and diagnostics
    <repo>/out/m5               MVMR instruments and estimates
    <repo>/out/g8               external replication output
    <repo>/out/sev              severity-stratified output
    <repo>/out/tables           manuscript tables built by 113_make_tables.py
    <repo>/out/figures          Fig 1 and Fig S1-S4 renders
    <repo>/figures              Fig 2 and Fig 3 renders

Large inputs (FinnGen R13, 1000G reference) are not committed.  Point the
environment variables at another location when they live outside the checkout:

    export ENDO_VTE_DATA=/scratch/endo_vte/data
    export ENDO_VTE_OUT=/scratch/endo_vte/out
    export ENDO_VTE_FIGURES=/scratch/endo_vte/figures
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

DATA_DIR = Path(os.environ.get("ENDO_VTE_DATA", REPO_ROOT / "data")).resolve()
OUT_DIR = Path(os.environ.get("ENDO_VTE_OUT", REPO_ROOT / "out")).resolve()
FIGURES_DIR = Path(os.environ.get("ENDO_VTE_FIGURES", REPO_ROOT / "figures")).resolve()
OUT_FIGURES_DIR = OUT_DIR / "figures"

# per-module output directories
LDSC_DIR = OUT_DIR / "ldsc"
MAGMA_DIR = OUT_DIR / "magma"
M3_DIR = OUT_DIR / "m3"
MR_DIR = OUT_DIR / "mr"
M5_DIR = OUT_DIR / "m5"
G8_DIR = OUT_DIR / "g8"
SEVERITY_DIR = OUT_DIR / "sev"
TABLE_DIR = OUT_DIR / "tables"

# primary FinnGen R13 inputs (GRCh38)
ENDO_SUMSTATS = DATA_DIR / "finngen_R13_N14_ENDOMETRIOSIS.gz"
VTE_SUMSTATS = DATA_DIR / "finngen_R13_I9_VTE.gz"

# LDSC reference files
LDSC_REF_LD = DATA_DIR / "ref_ld" / "LDscore."
LDSC_REF_W = DATA_DIR / "ref_w" / "1000G_Phase3_weights_hm3_no_MHC" / "weights.hm3_noMHC."
PLINK_1KG_DIR = DATA_DIR / "1000g_eur" / "1000G_EUR_Phase3_plink"

# external endometriosis meta-analysis (GWAS Catalog GCST90205183)
EXTERNAL_ENDO = DATA_DIR / "GCST90205183.h.tsv.gz"
