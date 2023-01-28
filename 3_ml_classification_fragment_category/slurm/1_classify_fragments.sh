#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_mea
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/3_ml_classification_fragment_category/src
#SBATCH -o /data5/deepro/starrseq/papers/results/3_ml_classification_fragment_category/slurm/logs/1_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/3_ml_classification_fragment_category/slurm/logs/1_err_%a.log
#SBATCH --exclude ramona,durga
#SBATCH --array 1-20
#SBATCH --spread-job

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

conda activate starrseq

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/3_ml_classification_fragment_category/slurm/files/1_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/3_ml_classification_fragment_category/src/1_classify_fragments.py $LINE

echo `date` ending job