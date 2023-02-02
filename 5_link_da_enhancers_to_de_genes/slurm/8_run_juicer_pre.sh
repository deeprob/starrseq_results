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

juicer_tools_path="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/Juicer/scripts/common/juicer_tools.jar"
stats_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/aligned/inter.txt"
graph_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/aligned/inter_hists.m"
restriction_sites="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/restriction_sites/hg38_MboI_mod.txt"
paired_contact_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/aligned/merged_nodups.txt"
out_hic_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/aligned/inter.hic"
chrom_sizes_file="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/Juicer/references/GRCh38.chrom.sizes.simple.sorted"
tmp_dir="./HIC_tmp"

################################################
#### setup environmental variables for java ####
################################################

mkdir -p ${tmp_dir}
HOME="/data5/deepro/tmp"
export HOME
# export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1""
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data5/deepro/tmp -Djava.util.prefs.systemRoot=/data5/deepro/tmp -Duser.home=/data5/deepro/tmp -Djava.library.path=/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common -Xmx128g -Xms49152m"
which java


# constructing hic matrix with KR normalization for 5000bp resolution as required by ABC input
java -Ddevelopment=false -Djava.awt.headless=true -jar ${juicer_tools_path} pre -s ${stats_file} -g ${graph_file} -f ${restriction_sites} -r 5000 -v ${paired_contact_file} -t ${tmp_dir} ${out_hic_file} ${chrom_sizes_file}

echo `date` ending job
