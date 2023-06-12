#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/src
#SBATCH -o /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/logs/3_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/logs/3_err_%a.log
#SBATCH --exclude laila,ramona
#SBATCH --array 1-46


source /opt/anaconda/bin/activate /data5/deepro/miniconda3/envs/starrseq

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/files/3_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/src/3_get_scores.py $LINE

echo `date` ending job
