# Code map

Script → output → manuscript location → value that must match after a rerun.
Section names refer to the v1.1 manuscript. All paths are relative to the repository root.

Verification rule: recompute the listed output from `data/`, then compare the anchor value with the manuscript; a mismatch means the repository snapshot and the manuscript are out of sync.

## 01 — sumstats download, munging and SNP tables

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `download_data.sh` | `data/finngen_R13_*.gz`, `1000G_Phase3_*.tgz` | Methods: Study design and data sources; Table I | Endometriosis 20,534 cases / 122,622 controls; VTE 27,928 / 472,258 |
| `pdl.py` | helper for `download_data.sh` | — | no manuscript value |
| `make_ldsc_sumstats.sh` | `data/endo_ldsc.sumstats.gz`, `data/vte_ldsc.sumstats.gz` | Methods: Genetic correlation analysis | LDSC N = 143,156 (endometriosis) and 500,186 (VTE) |

## 02 — LDSC main analysis (no archived script)

LDSC 1.0.1 was invoked directly. The logs are the only record.

| Output | Manuscript | Anchor value |
|---|---|---|
| `out/ldsc/h2_endo.log`, `h2_vte.log` | Results: Genetic correlation; Table I continued | h2 0.0699 (SE 0.0064) and 0.0185 (SE 0.0030); intercepts 1.139 and 1.121; attenuation ratios 0.372 and 0.354 |
| `out/ldsc/rg_endo_vte.log` | Results: Genetic correlation; Table I continued | rg 0.0443 (SE 0.0620), 95% CI −0.077 to 0.166, P = 0.475; cross-trait intercept 0.0086 (SE 0.0080) |

No `TableI_LDSC.tsv` is produced by the archived code; the LDSC block of Table I is transcribed from these logs.

## 03 — LDSC sensitivity (S5)

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `87_s5_ldsc_coag_excl.py` | `out/ldsc/s5_rg_exclusion_summary.tsv` | Results: Genetic correlation between endometriosis and VTE | ±1 Mb rg 0.048 (−0.065 to 0.162, P = 0.404); ±500 kb rg 0.051 (SE 0.0598, P = 0.390); gcov_int 0.0095 and 0.0105; change from baseline 0.007 |

## 04 — candidate regions and colocalization

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `80_m3_regions.py` | `out/m3/candidate_regions.tsv` | Results: A single high-confidence shared causal variant; Methods: Colocalization | 58 regions from 37 endometriosis and 46 VTE lead variants at ±500 kb; MHC kept separate |
| `81_m3_extract.py` | `out/m3/region_snps_endo.tsv.gz`, `region_snps_vte.tsv.gz` | Methods: Colocalization | region-level SNP extraction; 36 and 44 variants survive harmonisation |
| `82_m3_coloc.py` | `out/m3/coloc_results.tsv`, `coloc_sensitivity.tsv` | Results: A single high-confidence shared causal variant; Table IV | F5 PP.H4 0.994, PP.H3 0.006; PROCR PP.H4 0.502; ABO PP.H4 0.045, PP.H3 0.953; one region with PP.H4 ≥ 0.80 |

## 05 — p12 prior sensitivity

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `86_m3_p12_sens.py` | `out/m3/coloc_p12_sensitivity.tsv` | Results: A single high-confidence shared causal variant | F5 PP.H4 0.994 / 0.987 / 0.940 and PP.H3 0.006 / 0.013 / 0.060 at p12 = 1e-5 / 5e-6 / 1e-6; regions with PP.H4 ≥ 0.80: 1 / 1 / 1; PP.H3 ≥ 0.80: 7 / 7 / 8; PROCR PP.H4 0.502 / 0.335 / 0.092 |

## 06 — PLACO

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `83_m3_placo.py` | `out/m3/placo_all_scores.npz`, `placo_hits_p1e5.tsv.gz`, `placo_params.txt` | Results: A single high-confidence shared causal variant; Table IV (PLACO min P) | F5 minimum P 2.35 × 10⁻⁸²; ABO 2.15 × 10⁻⁵³ |
| `84_m3_placo_clump.py` | `out/m3/placo_loci.tsv` | Methods: Colocalization and pleiotropy analysis | clumping at r² < 0.1 within 1 Mb |
| `85_m3_summary.py` | `out/m3/m3_high_confidence.tsv`, `m3_region_summary.tsv`, `m3_report.txt` | Results: A single high-confidence shared causal variant; Table IV | one high-confidence shared region (F5) |

## 07 — MAGMA gene-set and permutation analysis

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `12_geneset_perm.py` | `out/geneset_enrichment.tsv`, `out/gene_stats.tsv.gz` | Results: Coagulation genes carry VTE-specific rather than shared signal; Fig 2B | VTE coagulation −log10 P = 13.8, BH-FDR q = 9.4 × 10⁻¹⁴; endometriosis P = 0.076; minimum-P permutation P = 0.0085 |
| `make_magma_pval.sh` | `data/endo_magma.pval.gz`, `vte_magma.pval.gz` | Methods: Gene-set and gene-level analysis | input preparation, no manuscript value |
| `make_snp_table.sh` | `data/snp_endo.tsv.gz`, `snp_vte.tsv.gz` | Methods: Gene-set and gene-level analysis | input preparation, no manuscript value |

## 08 — instruments and MR diagnostics

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `70_mr_extract.py` | `out/mr/endo_sig.tsv`, `vte_sig.tsv` | Methods: Instrument selection; Table II | 43 endometriosis and 84 VTE genome-wide significant loci; 36 and 44 clumped instruments |
| `71_mr_harmonize.py` | `out/mr/harmonised_forward.tsv`, `harmonised_reverse.tsv` | Table II; Table S2 | 36 and 44 instruments; mean F 60.9 and 123.1; minimum F 29.7 and 30.2 |
| `73_mr_diagnostics.py` | `out/mr/diagnostics_forward.tsv`, `diagnostics_reverse.tsv` | Results: No effect detected in the forward causal test; Fig S1 | forward I² = 0.35; per-SNP Wald ratios, Q contributions and MR-PRESSO outlier P |

## 09 — bidirectional MR

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `72_mr_analysis.py` | `out/mr/mr_results.tsv` | Results: The reverse direction suggests a protective association; Table III | reverse IVW OR 0.889 (0.840–0.940); weighted median 0.842 (0.748–0.947); unselected forward set (36 variants) OR 0.978 (0.768–1.245), I² = 0.98 |
| `74_mr_coag_sensitivity.py` | `out/mr/mr_coag_excluded.tsv` | Results: No effect detected in the forward causal test; The reverse direction; Table III | forward (34 variants) OR 1.017 (0.975–1.062), P = 0.432, I² = 0.35; MR-PRESSO OR 1.003 (0.968–1.039); reverse after coagulation exclusion (25 variants) IVW 0.887 (0.824–0.954), weighted median 0.798 (0.681–0.934) |

## 10 — multivariable MR

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `96_m5_build_union.py` | `out/m5/union_all.tsv`, `union_assoc.txt` | Methods: Multivariable Mendelian randomization | union of endometriosis, BMI and ever-smoking instruments |
| `97_m5_fetch.py` | `out/m5/instrument_matrix.tsv`, `instrument_effects_long.tsv` | Methods: Multivariable Mendelian randomization | aligned FinnGen R13 effects for each instrument |
| `98_m5_mvmr.py` | `out/m5/mvmr_results.tsv`, `mvmr_table.tsv` | Results: Forward multivariable models; Table III | BMI OR 1.47 (1.40–1.55); ever smoking 1.16 (1.09–1.24); conditional F 4.4 (BMI) and 19.3 (smoking); smoking-adjusted endometriosis OR 1.01 (0.973–1.05) |
| `99_m5_reverse_mvmr.py` | `out/m5/reverse_mvmr_results.tsv` | Results: Reverse-direction multivariable MR; Table S5 | univariable 0.936 (0.904–0.970, P = 2.4 × 10⁻⁴); +BMI 0.937 (0.881–0.997); +smoking 0.948 (0.892–1.007); +both 0.948 (0.891–1.009, P = 0.093); conditional F 71–111 (VTE), 3.4 (BMI), 1.5 (smoking); coagulation-excluded (27 variants) 0.915 (0.842–0.995) |

## 11 — F5 cross-phenotype lookup

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `60_fetch_fig3_data.py` | `out/fig3a_f5_spectrum.tsv`, `out/fig3b_coag_genes.tsv` | Results: Opposite effect directions across two phenotypic spectra; Fig 3A; Table S1A | DVT 3.28 (3.04–3.54); other embolism 2.75 (2.49–3.03); VTE 2.35 (2.23–2.47); PE 1.69 (1.57–1.82); ovarian cyst 0.92 (0.87–0.98); infertility 0.87 (0.80–0.94); endometriosis 0.79 (0.73–0.85); menorrhagia 0.69 (0.65–0.74) |

## 12 — sentinel genes, external replication and severity

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `62_coag_gene_sentinel.py` | `out/fig3b_sentinel_or.tsv` | Results: Coagulation genes carry VTE-specific rather than shared signal; Fig 3B; Table S1B | all nine genes stronger in VTE than endometriosis; four directionally opposite; FGG P = 0.061, PROCR P = 0.057; ABO endometriosis OR 0.957, VTE OR 0.788 |
| `90_g8_external_validation.py` | `out/g8/g8_f5_lookup.tsv`, `g8_reverse_mr.tsv`, `g8_report.txt` | Results: The F5 association is directionally concordant in the external data; The reverse direction; Table V | external F5 OR 0.866 (0.798–0.939, P = 5.18 × 10⁻⁴); proxy 0.853 (0.788–0.924, P = 1.03 × 10⁻⁴); reverse external (38 variants) OR 0.926 (0.880–0.974), weighted median 0.889 (0.782–1.011, P = 0.073) |
| `99_severity_pe_dvt.py` | `out/sev/severity_instruments.tsv`, `severity_mr.tsv`, `severity_specific.tsv`, `severity_mvmr.tsv`, `severity_vte_split.tsv` | Results: Severity stratification does not support diagnostic competition; Table VI | PE-instrumented (27) OR 0.953 (0.876–1.037); DVT-instrumented (46) 0.944 (0.900–0.991); DVT-specific (34) 0.865 (0.821–0.910); PE-specific (6) 0.914 (0.758–1.103); PE/DVT ratio 1.06 (0.87–1.28, P = 0.573) |
| `100_severity_coag_sens.py` | stdout (MR-PRESSO global P; severity estimates after coagulation-window removal) | Results: No effect detected in the forward causal test; Severity stratification; Table VI | forward coagulation-excluded MR-PRESSO P = 0.880, I² = 0 |
| `114_severity_wmedian.py` | `out/sev/severity_wmedian_extra.tsv` | Results: Severity stratification; Table VI | DVT weighted median 0.896 (0.811–1.017); DVT-specific weighted median 0.854 (0.781–0.976) |

## 13 — endpoint table and figures

| Script | Output | Manuscript | Anchor value |
|---|---|---|---|
| `20_make_fig2.py` | `figures/fig2_rg_geneset.png|pdf|tiff` | Fig 2A, 2B; Results: Genetic correlation; Coagulation genes | rg 0.0443; VTE coagulation −log10 P = 13.8 |
| `21_make_fig3.py` | `figures/fig3_f5_coag.png|pdf|tiff` | Fig 3A, 3B; Results: Opposite effect directions; Coagulation genes | Fig 3A ORs as in module 11; Fig 3B sentinel ORs as in module 12 |
| `111_figS1_S2.py` | `out/figures/FigS1_diagnostic.*`, `FigS2_MR_full.*` | Fig S1, Fig S2 | forward and reverse MR estimates and diagnostics |
| `112_figS4_locuszoom.py` | `out/figures/FigS4_LocusZoom.*` | Fig S4 | F5 and ABO regional LD (1000G EUR) |
| `113_make_tables.py` | `out/tables/TableI_data_sources.tsv`, `TableII_instruments.tsv`, `TableIV_colocalization.tsv`, `TableV_external_validation.tsv`, `TableVI_severity.tsv` | Tables I, II, IV, V, VI | Table I: h2 0.0699 / 0.0185, rg 0.0443; Table II: 36 / 44 instruments; Table IV: F5 PP.H4 0.994; Table V: external F5 0.866; Table VI: PE 0.953, DVT 0.944 |
| `115_fig1_workflow.py` | `out/figures/Fig1_study_design.*` | Fig 1 | study design and workflow |
| `116_fig1_methods.py` | `out/figures/Fig1_study_design.*` (methods-only variant) | Fig 1 (methods version) | no numeric labels; the script asserts that no digit is displayed |

## Not covered by an archived script

| Item | Where it appears | Recorded value |
|---|---|---|
| LDSC h2/rg driver (module 02) | Table I continued; Results: Genetic correlation | logs in `out/ldsc/` |
| liftOver and QC (module 01) | Methods: Quality control and harmonisation | munged LDSC files are the starting point |
| Three prevalence tiers, MHC-excluded LDSC, split-half (seed 20261108) | Methods: Sensitivity analyses | not archived |
| sdY prior sensitivity (0.10–0.30) | Results: A single high-confidence shared causal variant | PP.H4 0.993–0.994 |
| Risteys/PheWeb endpoint definitions | Table S3 | not archived |
| Fig S3 (PLACO panel) | Results: A single high-confidence shared causal variant | no producer script |
| `TableI_LDSC.tsv` | checklist item for module 02 | never generated; LDSC block transcribed from logs |
