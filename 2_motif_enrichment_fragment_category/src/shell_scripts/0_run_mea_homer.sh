#!/bin/bash
set -ue

region_file=$1
reference_genome=$2
background_region=$3
output_file=$4
threads=$5

findMotifsGenome.pl $region_file $reference_genome $output_file -bg $background_region -p $threads
