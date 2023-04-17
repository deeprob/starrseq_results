#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_abc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=200G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/11_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/11_err_%a.log
#SBATCH --exclude ramona,durga
#SBATCH --array 3-4
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

conda activate final-abc-env

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/files/11_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/ABC-Enhancer-Gene-Prediction/src/predict.py $LINE

echo `date` ending job