#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_act
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/src
#SBATCH -o /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/logs/1_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/logs/1_err.log
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

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/slurm/files/1_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/src/1_store_categorized_fragments.py $LINE

echo `date` ending job
