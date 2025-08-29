#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_act
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/src
#SBATCH -o /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/logs/2_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/logs/2_err.log
#SBATCH --exclude ramona,durga,qingyu
#SBATCH --array 1-6


source /opt/anaconda/bin/activate /data5/deepro/miniconda3/envs/starrseq

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/files/2_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/src/2_create_direct_binding_fragment.py $LINE

echo `date` ending job
