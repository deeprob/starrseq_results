#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_target
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/src
#SBATCH -o /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/slurm/logs/14_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/slurm/logs/14_err.log
#SBATCH --exclude ramona,durga


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

export MPLCONFIGDIR="/data5/deepro/tmp"

echo `date` starting job on $HOSTNAME

python /data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/src/14_create_da_de_table.py

echo `date` ending job
