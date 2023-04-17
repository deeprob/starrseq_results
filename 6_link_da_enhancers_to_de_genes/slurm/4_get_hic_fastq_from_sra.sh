#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_hic
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes
#SBATCH -o /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/4_out_%a.log
#SBATCH -e /data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/slurm/logs/4_err_%a.log
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

conda activate ncbi_download

echo `date` starting job on $HOSTNAME

# download sra accession from the project page link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP018571&o=acc_s%3Aa
# manually remove accession numbers which are not relevant - only keep Hi-C, HEK293 siRNA Control SRA runs
# esearch -db sra -query SRP018571 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > SRR_Acc_List.txt

output_dir="/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t"
# prefetch files
prefetch --option-file ${output_dir}/SRR_Acc_List.txt --output-directory $output_dir

# manually edit the accession list to select only required datsets

# download fastq files for SRR data
cat ${output_dir}/SRR_Acc_List.txt | xargs -I {} fasterq-dump --split-files -t /data5/deepro/tmp ${output_dir}/{}/{}.sra -O $output_dir -e 8

fastq_dir="${output_dir}/fastq"
mv ${output_dir}/*.fastq fastq_dir

echo `date` ending job
