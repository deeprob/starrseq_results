#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_juice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=74
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=32G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/6_out.log 
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/6_err.log
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

conda activate hic_pre

echo `date` starting job on $HOSTNAME

juicer_dir="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6"
top_dir="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t"
reference_genome="/data5/deepro/genomes/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
restriction_sites="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/restriction_sites/hg38_MboI_mod.txt"
chrom_sizes="/data5/deepro/genomes/hg38/GRCh38.chrom.sizes.simple.sorted"

################################################
#### setup environmental variables for java ####
################################################

HOME="/data5/deepro/tmp"
export HOME
# export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1""
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data5/deepro/tmp -Djava.util.prefs.systemRoot=/data5/deepro/tmp -Duser.home=/data5/deepro/tmp -Djava.library.path=/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common -Xmx32g -Xms4096m"
which java

${juicer_dir}/scripts/juicer.sh -d ${top_dir} -D ${juicer_dir} -z ${reference_genome} -p ${chrom_sizes} -y ${restriction_sites} -t 64 -g GRCh38 -S early

echo `date` ending job
