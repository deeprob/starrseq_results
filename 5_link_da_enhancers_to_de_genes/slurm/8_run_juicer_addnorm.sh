#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_juice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=200G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/8_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/8_err.log
#SBATCH --nodelist sarah

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

juicer_tools_path="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/scripts/common/juicer_tools.jar"
out_hic_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/aligned/inter.hic"

################################################
#### setup environmental variables for java ####
################################################

HOME="/data5/deepro/tmp"
export HOME
# export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1""
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data5/deepro/tmp -Djava.util.prefs.systemRoot=/data5/deepro/tmp -Duser.home=/data5/deepro/tmp -Djava.library.path=/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common -Xmx128g -Xms49152m"
which java


# constructing hic matrix with KR normalization as required by ABC input
java -Ddevelopment=false -Djava.awt.headless=true -jar ${juicer_tools_path} addNorm -F -w 5000 ${out_hic_file}

echo `date` ending job
