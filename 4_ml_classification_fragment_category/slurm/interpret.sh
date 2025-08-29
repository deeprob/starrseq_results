#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=cc_interpret 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data7/deepro/starrseq/4_ml_classification_fragment_category # TODO: set dir to data dir
#SBATCH -o /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/out_interpret_%a.log # TODO: set slurm output file
#SBATCH -e /data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/logs/err_interpret_%a.log # TODO: set slurm input file
#SBATCH --nodelist=sarah
#SBATCH --array 1-8

export HOME="/data7/deepro/starrseq/4_ml_classification_fragment_category/data/tmp"
echo `date` starting job on $HOSTNAME
source /opt/anaconda/bin/activate /data7/deepro/miniconda3/envs/dlinterpret

motif_path="/data7/deepro/starrseq/4_ml_classification_fragment_category/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

smap_file="/data7/deepro/starrseq/4_ml_classification_fragment_category/slurm/files/interpret.txt"

ohe_path=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" $smap_file)
attr_path=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" $smap_file)
out_path=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" $smap_file)
result_dir=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$4}" $smap_file)

modisco motifs -s $ohe_path -a $attr_path -n 2000 -o $out_path -w 500
modisco report -i $out_path -o $result_dir -s $result_dir -t -n 3 -m $motif_path
modisco meme -i $out_path -o ${result_dir} -s $result_dir -t -n 3 -m $motif_path

echo `date` ending job on $HOSTNAME
