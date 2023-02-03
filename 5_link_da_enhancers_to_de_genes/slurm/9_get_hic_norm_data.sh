#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_hic
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=128G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/8_out_%a.log 
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/8_err_%a.log 
#SBATCH --nodelist sarah
#SBATCH --array 1

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

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/files/8_smap.txt)

juicer_tools_path="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/scripts/common/juicer_tools.jar"

################################################
#### setup environmental variables for java ####
################################################

HOME="/data5/deepro/tmp"
export HOME
# export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1""
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data5/deepro/tmp -Djava.util.prefs.systemRoot=/data5/deepro/tmp -Duser.home=/data5/deepro/tmp -Djava.library.path=/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common -Xmx128g -Xms49152m"

echo $LINE
python /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/ABC-Enhancer-Gene-Prediction/src/juicebox_dump.py --juicebox "/data5/deepro/miniconda3/envs/hic_pre/bin/java -Ddevelopment=false -Djava.awt.headless=true -jar ${juicer_tools_path}" $LINE

echo `date` ending job
