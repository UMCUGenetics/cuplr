#!/bin/bash
#SBATCH --output=train_final_model.o
#SBATCH --job-name=train_final
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=15

wd=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/
mkdir -p $wd/final/

if [[ ! -f $wd/final/job.done ]]; then
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript $wd/do_train.R $wd/features/features.txt.gz '' $wd/final/model.rds 1 && touch $wd/final/job.done
EOF
else
echo Skipping. Done file exists: $wd/final/job.done
fi

