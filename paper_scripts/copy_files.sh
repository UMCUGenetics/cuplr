#!/bin/bash

paper_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/paper_scripts/
models_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/

## Main model ================================
main_model_in_dir=$models_dir/0.21a_DR104-update5_pcawgSoftFilt2/
main_model_out_dir=$paper_dir/01_main_model/
mkdir -p $main_model_out_dir

cp $main_model_in_dir/do_train.R $main_model_out_dir/00_rf_train_template.R
cp $main_model_in_dir/setup_training.R $main_model_out_dir/01_setup_training_data.R
cp $main_model_in_dir/spawn_cv_jobs.R $main_model_out_dir/02_spawn_cv_jobs.R
cp $main_model_in_dir/train_final_model.sh $main_model_out_dir/03_train_final_model.sh
cp $main_model_in_dir/training_report.R $main_model_out_dir/04_analysis.R

## Feature exclusion model ================================
feat_excl_model_in_dir=$models_dir/0.21c_featGroupExcl/
feat_excl_model_out_dir=$paper_dir/02_feature_exclusion_models/
mkdir -p $feat_excl_model_out_dir

cp $feat_excl_model_in_dir/do_train.R $feat_excl_model_out_dir/00_rf_train_template.R
cp $feat_excl_model_in_dir/setup_training.R $feat_excl_model_out_dir/01_setup_training_data.R
cp $feat_excl_model_in_dir/setup_excl_feat_cv.R $feat_excl_model_out_dir/02_setup_excl_feat_cv.R
cp $feat_excl_model_in_dir/excl_report.R $feat_excl_model_out_dir/03_compare_feature_exclusion_models.R

## Clonal muts model ================================
clonal_muts_model_in_dir=$models_dir/0.21b_DR104-update5_pcawgSoftFilt2_clonal/
clonal_muts_model_out_dir=$paper_dir/03_clonal_mutations_model/
mkdir -p $clonal_muts_model_out_dir

cp $clonal_muts_model_in_dir/do_train.R $clonal_muts_model_out_dir/00_rf_train_template.R
cp $clonal_muts_model_in_dir/setup_training.R $clonal_muts_model_out_dir/01_setup_training_data.R
cp $clonal_muts_model_in_dir/spawn_cv_jobs.R $clonal_muts_model_out_dir/02_spawn_cv_jobs.R
cp $clonal_muts_model_in_dir/train_final_model.sh $clonal_muts_model_out_dir/03_train_final_model.sh
cp $clonal_muts_model_in_dir/training_report.R $clonal_muts_model_out_dir/04_training_report.R
cp $clonal_muts_model_in_dir/compare_all_vs_clonal_cuplr.R $clonal_muts_model_out_dir/05_compare_main_vs_clonal_only_model.R

## Hartwig or PCAWG only model ================================
hartwig_pcawg_only_model_in_dir=$models_dir/0.21d_HMF_or_PCAWG_train/
hartwig_pcawg_only_model_out_dir=$paper_dir/04_hartwig_or_pcawg_only_model/
mkdir -p $hartwig_pcawg_only_model_out_dir

cp $hartwig_pcawg_only_model_in_dir/do_train.R $hartwig_pcawg_only_model_out_dir/00_rf_train_template.R
cp $hartwig_pcawg_only_model_in_dir/setup_training.R $hartwig_pcawg_only_model_out_dir/01_setup_training_data.R
cp $hartwig_pcawg_only_model_in_dir/spawn_cv_jobs.R $hartwig_pcawg_only_model_out_dir/02_spawn_cv_jobs.R
cp $hartwig_pcawg_only_model_in_dir/compare_models.R $hartwig_pcawg_only_model_out_dir/03_compare_hartwig_vs_pcawg_models.R



