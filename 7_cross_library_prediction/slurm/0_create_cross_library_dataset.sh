#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/7_cross_library_prediction/src
#SBATCH -o /data5/deepro/starrseq/papers/results/7_cross_library_prediction/slurm/logs/0_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/7_cross_library_prediction/slurm/logs/0_err_%a.log
#SBATCH --exclude laila,ramona
#SBATCH --array 1


source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/dlenv

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/7_cross_library_prediction/slurm/files/0_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/7_cross_library_prediction/src/0_create_cross_library_dataset.py $LINE

echo `date` ending job
