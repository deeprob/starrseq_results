#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=cc_train 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data7/deepro/starrseq/4_ml_classification_fragment_category # TODO: set dir to data dir
#SBATCH -o /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/out_eval_%a.log # TODO: set slurm output file
#SBATCH -e /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/err_eval_%a.log # TODO: set slurm input file
#SBATCH --nodelist=laila
#SBATCH --gpus=1
#SBATCH --array 1-8

export HOME="/data7/deepro/starrseq/4_ml_classification_fragment_category/data/tmp"
echo `date` starting job on $HOSTNAME
source /opt/anaconda/bin/activate /data7/deepro/miniconda3/envs/biopeaker

peaker_path="/data7/deepro/pipelines/biopeaker/src/peaker.py"
scores_path="/data7/deepro/starrseq/4_ml_classification_fragment_category/src/1_get_scores.py"

smap_file="/data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/files/evaluate.txt"

dataset_path=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" $smap_file)
save_dir=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" $smap_file)
# task_name=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" $smap_file)

genome_fasta="/data5/deepro/genomes/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

pred_file="${save_dir}resnet_mlp.csv.gz"
score_file="${save_dir}scores.csv"


python $peaker_path $dataset_path $genome_fasta $save_dir --test --interpreter ig # --label_name $task_name
python $scores_path $pred_file $score_file

echo `date` ending job on $HOSTNAME
