#!/bin/bash
# set -ue # don't set ue because it interrupts conda activate


# activate conda environment
source /opt/anaconda/bin/activate /data5/deepro/miniconda3/envs/gseanew


gene_file=$1
gseaout_file=$2
keggout_file=$3
gseafigout_file=$4
keggfigout_file=$5
tmp_dir=$6

Rscript /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/src/scripts/enrich.R $gene_file $gseaout_file $keggout_file $gseafigout_file $keggfigout_file $tmp_dir
