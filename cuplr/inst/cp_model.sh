#!/bin/bash

out_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/inst/
model_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/

cp $model_dir/final/model.rds $out_dir/model.rds
cp $model_dir/report/prob_calib_curves.txt $out_dir/prob_calib_curves.txt
