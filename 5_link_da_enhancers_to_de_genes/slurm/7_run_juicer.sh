#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_juice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=30G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/7_out.log
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes//slurm/logs/7_err.log
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

juicer_dir="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/Juicer"
top_dir="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t"

${juicer_dir}/scripts/juicer.sh -d ${top_dir} -D ${juicer_dir} -z ${juicer_dir}/references/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -p ${juicer_dir}/references/GRCh38.chrom.sizes.simple.sorted -y ${juicer_dir}/restriction_sites/hg38_MboI.txt -t 32 -T 32 -g GRCh38

echo `date` ending job
