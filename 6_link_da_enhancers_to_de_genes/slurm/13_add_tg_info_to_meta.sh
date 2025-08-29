#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_target
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/src
#SBATCH -o /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/slurm/logs/13_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/slurm/logs/13_err.log
#SBATCH --exclude ramona,durga


export HOME="/data5/deepro/tmp"

source /opt/anaconda/bin/activate /data5/deepro/miniconda3/envs/starrseq

echo `date` starting job on $HOSTNAME

python /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/src/13_add_tg_info_to_meta.py

echo `date` ending job
