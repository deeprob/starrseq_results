#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rnaseq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/4_compare_expression_ko_vs_wt
#SBATCH -o /data5/deepro/starrseq/papers/results/4_compare_expression_ko_vs_wt/slurm/logs/out_r2c_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/4_compare_expression_ko_vs_wt/slurm/logs/err_r2c_%a.log
#SBATCH --exclude=ramona
#SBATCH --array 1-8

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


conda activate rnaseq

##### This script runs rnaseq reads to counts pipeline on prabhu lab data ####

echo `date` starting job on $HOSTNAME
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/4_compare_expression_ko_vs_wt/slurm/files/0_smap.txt)

python /data6/deepro/computational_pipelines/rnaseq_analysis/src/0_reads_to_counts.py $LINE

echo `date` ending job on $HOSTNAME
