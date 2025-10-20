# Quantifying the impact of mutations on enhancer dynamics
Codebase for manuscript titled: "Quantifying the impact of mutations on enhancer dynamics" by Das and Banerjee et al.

# Note
Raw data processing and quality control pipeline of the processed data analyzed in this repo is described here: https://github.com/deeprob/starrseq_reproducibility

# Repository guidelines
This repo contains step-by-step guide to the analysis carried out for " ". A description of each sub dir is as follows:

- **0_quality_control**: STARR-seq data quality control scripts.

- **1_compare_activity_ko_vs_wt**: Differential enhancer activity calling scripts.

- **2_categorize_fragment_on_activity**: Defining categories of differentially active enhancers.

- **3_motif_enrichment_fragment_category**: Motif enrichment analysis of various categories of enhancers.

- **4_ml_classification_fragment_category**: Deep Learning models to classify fragments based on their enhancer activity.

- **5_compare_expression_ko_vs_wt**: Differential gene expression scripts.

- **6_link_da_enhancers_to_de_genes**: Mapping enhancers to target genes scripts.

