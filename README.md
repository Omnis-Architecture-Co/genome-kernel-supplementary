# Genome Kernel Supplementary Materials

**"The human genome satisfies the formal properties of a deterministic computational kernel"**

Jasmine Levy • OMNIS Architecture Co.

---

## Overview

This repository contains all source data, validation scripts, and analysis code for the manuscript. Every statistical claim in the paper can be reproduced from these materials.

## Quick Start

```bash
# Install dependencies
pip install numpy scipy pandas

# Run the complete validation suite (21 tests, ~35 seconds)
cd source_data/
python ../scripts/validation_suite.py

# Run just the full-reproduction tests
python ../scripts/validation_suite.py --full

# Run a specific test
python ../scripts/validation_suite.py VAL-03

# List all available tests
python ../scripts/validation_suite.py --list
```

## Repository Structure

```
genome-kernel-supplementary/
│
├── README.md
│
├── scripts/
│   ├── validation_suite.py          # All 21 statistical tests in one script
│   ├── upstream_pic_analysis.py     # 100-gene PIC carrier protein analysis
│   ├── permutation_test.py          # Multi-seed permutation test
│   ├── worked_examples.py           # TP53 and OR7D4 worked examples
│   ├── fetch_genes.py               # Gene sequence retrieval from Ensembl
│   └── fetch_coords_and_controls.py # Random control generation
│
├── source_data/
│   ├── vocabulary_human_1932words.csv
│   ├── vocabulary_mouse_1117words.csv
│   ├── vocabulary_zebrafish_2403words.csv
│   ├── vocabulary_fly_937words.csv
│   ├── vocabulary_yeast_123words.csv
│   ├── vocabulary_ecoli_7words.csv
│   ├── vocabulary_celegans_430words.csv
│   ├── vocabulary_arabidopsis_534words.csv
│   ├── vocabulary_halobacterium_3words.csv
│   ├── programs_annotated_4936.csv
│   ├── dispatch_matrix_cross_chromosome.csv
│   ├── chromosome_roles.csv
│   ├── gene_departments.csv
│   ├── execution_trace_summary.csv
│   ├── execution_trace_hop1.csv
│   ├── primitive_annotations_116.csv
│   ├── TableS5_codon_encoding.csv
│   ├── TableS_department_classification.csv
│   ├── TableS_input_specifications.csv
│   ├── 100gene_PIC_analysis.csv
│   ├── gene_manifest.csv
│   └── gene_panel_coordinates.csv
│
├── validation_results/
│   ├── VAL-ENC-001_encoding_null_model.json
│   ├── VAL-CON-001_convergence_null_model.json
│   ├── VAL-PRM-001_primitive_recurrence.json
│   ├── VAL-NET-001_dispatch_hub_null.json
│   ├── VAL-XSP-001_cross_species_tau.json
│   ├── VAL-PEEL-ADDENDUM_results.json
│   ├── issue2_encoding_permutations.json
│   ├── isoform_enc001_canonical.json
│   ├── issue6_isoform_collapsed.json
│   ├── issue8_parameter_sensitivity.json
│   ├── chr19_density_normalization.json
│   ├── depmap_essentiality_results.json
│   ├── dispatch_vs_ppi_results.json
│   ├── full_pipeline_negative_controls.json
│   ├── issue17_confidence_robustness.json
│   ├── issue3_exonic_convergence.json
│   ├── vocab_vs_blast_results.json
│   ├── vocab_vs_blast_essentiality_results.json
│   ├── vocab_residual_essentiality_results.json
│   ├── pfam_comparison_results.json
│   ├── shuffled_genome_control_results.json
│   └── chrm_independence_results.json
│
├── supplementary_tables/
│   ├── Table_S1_permutation_invariance.csv
│   ├── Table_S2_isoform_collapse.csv
│   ├── Table_S3_parameter_sensitivity.csv
│   ├── Table_S4_gene_density_normalization.csv
│   ├── Table_S5_depmap_essentiality.csv
│   ├── Table_S6_string_validation.csv
│   ├── Table_S7_sequence_order.csv
│   ├── Table_S8_negative_control.csv
│   ├── Table_S9_token_sensitivity.csv
│   ├── Table_S10_priority_robustness.csv
│   └── Table_S11_residual_essentiality.csv
│
├── supplementary_figures/
│   ├── Figure5_upstream_token_analysis.png
│   ├── VALENC001_byte_distribution.png
│   ├── VALPRM001_recurrence_distribution.png
│   ├── figure3_cross_species_9sp.png
│   ├── VALPEELADDENDUM_comparison.png
│   └── VALDICT001_v6c_layered_peel_summary.png
│
├── supplementary_text/
│   ├── Supplementary_Methods.txt
│   ├── Supplementary_Note1_Worked_Examples.txt
│   └── README.txt
│
└── File_Manifest.csv
```

## Validation Suite

The `validation_suite.py` script reproduces all 21 statistical tests from the paper:

| ID | Test | Mode | Key Result |
|---|---|---|---|
| VAL-01 | Encoding null model | Verify | t=5.88, p=4.8×10⁻⁹ |
| VAL-02 | Convergence null model | Verify | z=86.4, p<0.001 |
| VAL-03 | Primitive recurrence | **Full** | z=52.4, p<0.002 |
| VAL-04 | Dispatch hub structure | **Full** | Gini z=18.1, p<0.001 |
| VAL-05 | Cross-species conservation | **Full** | τ=-0.449, p=0.0001 |
| ROB-01 | Encoding permutations | Verify | CV=2.0% |
| ROB-02 | Isoform collapse | Verify | Both significant |
| ROB-03 | Parameter sensitivity | Verify | CV=0.39 |
| ROB-04 | Gene density normalization | **Full** | chr19 rank 1 |
| ROB-05 | DepMap essentiality | **Full** | η²=0.190 |
| ROB-06 | STRING PPI validation | **Full** | 5.0× enrichment |
| ROB-07 | Sequence order dependence | Verify | 6,893× at 4 bytes |
| ROB-08 | Negative controls | Verify | 0 programs, 0 edges |
| ROB-09 | Token sensitivity | Verify | 12 configs stable |
| ROB-10 | Priority robustness | Verify | ρ=0.988 |
| ROB-11 | Residual essentiality | Verify | η²=0.024 |
| ROB-12 | Progressive peel | Verify | All sig at L2 |
| ROB-13 | Vocab vs BLAST | Verify | 11.9% accuracy |
| ROB-14 | Pfam comparison | Verify | Sub-domain resolution |
| ROB-15 | Exonic convergence | Verify | 68.5% overlap |
| V2-04 | Fisher's combined test | **Full** | p=5.36×10⁻¹⁶ |

**Full** = re-runs computation from source data. **Verify** = confirms published statistics from result JSONs.

## External Data

Two external datasets are required for full reproduction of ROB-05 and ROB-06:

- **DepMap 25Q3**: Download `CRISPRGeneEffect.csv` from [depmap.org](https://depmap.org/portal/data_page/?tab=allData&release=DepMap+Public+25Q3). Place in `source_data/external_data/`.
- **STRING v12.0**: Download `9606.protein.links.v12.0.txt.gz` from [string-db.org](https://string-db.org/cgi/download?species_text=Homo+sapiens). Decompress and place in `source_data/external_data/`.

These datasets were not used during vocabulary construction and serve as independent biological validations.

## Reproducibility

All random number generators use seed 42. The pipeline is deterministic: re-execution on identical input produces byte-identical output. MD5 checksums for all input files are provided in `File_Manifest.csv`.

## Requirements

- Python 3.9+
- NumPy, SciPy, Pandas

## Citation

If you use this data or code, please cite:

> Levy, J. (2026). The human genome satisfies the formal properties of a deterministic computational kernel. *Nature* [submitted].

## License

Source data and analysis code are provided for research purposes under CC BY 4.0. The OMNIS encoding pipeline is subject to a separate provisional patent (see manuscript Competing Interests).
