#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=peaks 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data7/deepro/starrseq/4_ml_classification_fragment_category # TODO: set dir to data dir
#SBATCH -o /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/out_peak_%a.log # TODO: set slurm output file
#SBATCH -e /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/err_peak_%a.log # TODO: set slurm input file
#SBATCH --nodelist=laila
#SBATCH --array 1-8

export HOME="/data7/deepro/starrseq/4_ml_classification_fragment_category/data/tmp"
echo `date` starting job on $HOSTNAME
source /opt/anaconda/bin/activate /data7/deepro/miniconda3/envs/biopeaker

peak_path="/data7/deepro/starrseq/4_ml_classification_fragment_category/src/0_generate_peak_files.py"
dataset_create_path="/data7/deepro/pipelines/biopeaker/src/create_dataset.py"

smap_file="/data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/files/peaks.txt"

peak_file=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" $smap_file)
master_file=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" $smap_file)
save_dir=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" $smap_file)

echo $peak_file
echo $master_file
echo $save_dir

peak_bed_path="${save_dir}peaks.bed.gz"
nonpeak_bed_path="${save_dir}not_peaks.bed.gz"
dataset_path="${save_dir}resnet_1000.h5"

python $peak_path $peak_file $master_file $save_dir
python $dataset_create_path $peak_bed_path $nonpeak_bed_path $dataset_path --split_ratio 70 15 15

echo `date` ending job on $HOSTNAME
