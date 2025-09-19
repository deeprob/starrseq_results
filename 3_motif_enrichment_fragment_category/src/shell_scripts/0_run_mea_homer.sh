#!/bin/bash
set -ue

region_file=$1
reference_genome=$2
background_region=$3
output_dir=$4
threads=$5
motif_file=$6

findMotifsGenome.pl $region_file $reference_genome $output_dir -bg $background_region -mknown $motif_file -p $threads -size given -nomotif
