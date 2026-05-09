SUPPLEMENTARY INFORMATION PACKAGE
==================================================

OMNIS Architecture: A Computational Vocabulary Encoded in the Human Proteome

Contents:

Supplementary_Methods.txt
  Fifteen individually numbered robustness analyses (i-xv) plus S11 residual analysis.

Supplementary_Note1_Worked_Examples.txt
  Complete encoding pipeline walkthroughs for Nucleolin (protein pathway) and
  INS locus (DNA pathway), referenced in Methods.

Supplementary_Tables/
  Table_S1_permutation_invariance.csv     - 6 encoding permutations, vocab sizes
  Table_S2_isoform_collapse.csv           - ENC-001 + CON-001 canonical-only tests
  Table_S3_parameter_sensitivity.csv      - Window, threshold, enrichment configs
  Table_S4_gene_density_normalization.csv  - 24 chromosomes, gene-density-corrected ratios
  Table_S5_depmap_essentiality.csv         - Department essentiality (DepMap 25Q3)
  Table_S6_string_validation.csv           - STRING PPI enrichment (3 groups)
  Table_S7_sequence_order.csv              - Shuffled genome, word-length stratified
  Table_S8_negative_control.csv            - Full pipeline on 4 control conditions
  Table_S9_token_sensitivity.csv           - Token-level sensitivity analysis
  Table_S10_priority_robustness.csv        - Priority weighting configs
  Table_S11_residual_essentiality.csv      - Residual essentiality beyond BLAST-proxy

Supplementary_Figures/
  Figures referenced in the main text and supplementary methods.

Source_Data/
  All nine species vocabulary CSVs and supporting data CSVs:
  Human, Mouse, Zebrafish, Fly, Yeast, E. coli,
  C. elegans, Arabidopsis, H. salinarum

Validation_JSONs/
  Machine-readable validation outputs from all robustness analyses.
  These are the raw data from which supplementary tables were generated.

Table Numbering:
  Tables S1-S11 are contiguous (no gaps). S11 was previously numbered S16.