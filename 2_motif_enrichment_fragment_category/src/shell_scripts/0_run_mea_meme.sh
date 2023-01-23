#!/bin/bash
set -ue

region_fasta_file=$1
background_region_fasta_file=$2
motif_file=$3
output_folder=$4

ame --oc $output_folder --control $background_region_fasta_file $region_fasta_file $motif_file
