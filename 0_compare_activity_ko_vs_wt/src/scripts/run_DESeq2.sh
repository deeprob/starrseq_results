#!/bin/bash
# set -ue # don't set ue because it interrrupts with conda activate

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

# activate conda environment
conda activate starrseq-da


infile=$1
designmat=$2
outfile=$3
echo $PWD

Rscript /data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/src/scripts/deseq2.R $infile $designmat $outfile
