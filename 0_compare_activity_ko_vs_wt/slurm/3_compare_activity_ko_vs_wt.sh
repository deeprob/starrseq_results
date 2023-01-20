#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_da_peaks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/src
#SBATCH -o /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/slurm/logs/3_out_0_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/slurm/logs/3_err_0_%a.log
#SBATCH --exclude durga,ramona
#SBATCH --array 1-7

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
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/slurm/files/3_smap.txt)

echo $LINE
python /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/src/3_compare_activity_ko_vs_wt.py $LINE

echo `date` ending job
