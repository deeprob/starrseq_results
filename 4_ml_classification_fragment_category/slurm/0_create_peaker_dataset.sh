#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/src
#SBATCH -o /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/logs/0_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/logs/0_err_%a.log
#SBATCH --nodelist laila
#SBATCH --array 1-23


source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/dlenv

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/4_ml_classification_fragment_category/slurm/files/0_smap.txt)

echo $LINE
python /data6/deepro/computational_pipelines/biopeaker/src/create_dataset.py $LINE

echo `date` ending job
