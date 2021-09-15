#!/bin/bash
#SBATCH --job-name=train_final
#SBATCH --output=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/cuplr/models/0.18a_newLabels/final/slurm.out
#SBATCH --time=3:00:00
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=15

wd=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/cuplr/models/0.18a_newLabels/
mkdir -p $wd/final/

if [[ ! -f $wd/final/job.done ]]; then
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript $wd/do_train.R $wd/features/features.rds '' $wd/final/model.rds 4 && touch $wd/final/job.done
EOF
else
echo Skipping. Done file exists: $wd/final/job.done
fi

