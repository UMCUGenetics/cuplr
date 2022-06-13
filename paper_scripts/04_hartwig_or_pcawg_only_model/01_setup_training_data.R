options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')

#devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/commonUtils/'))
#devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))

model_dir <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21d_HMF_or_PCAWG_train/')
out_dir <- paste0(model_dir,'/features/'); dir.create(out_dir, recursive=T, showWarnings=F)

## Load metadata ================================
file.copy(
   paste0(base_dir,'/CUPs_classifier/processed/metadata/12_20211210/04_metadata_training.txt'),
   paste0(out_dir,'/metadata_training.txt'),
   overwrite=T
)

metadata <- read.delim(paste0(out_dir,'/metadata_training.txt'))

## Helper functions ================================
getMetadata <- function(sample.names, keys='cancer_type'){
   metadata[match(sample.names, metadata$sample_id), keys]
}

## Prepare features ================================
message('Preparing features...')
df_features <- read.delim(
   paste0(base_dir,'/CUPs_classifier/processed/features/_all/09_DR104-update5_pcawgSoftFilt2/all_features.txt.gz'),
   check.names=F
)

## Add cancer type
df_features <- cbind(
   response=getMetadata(rownames(df_features), 'cancer_type'),
   df_features
)
df_features$response <- as.factor(df_features$response)

## Select training samples ================================
## Since held out samples will either be HMF only or PCAWG only, we can use both training and holdout samples for training
sample_whitelist <- metadata[metadata$is_training_sample | metadata$is_holdout_sample,'sample_id']

training_samples <- rownames(df_features)[rownames(df_features) %in% sample_whitelist]
training_samples <- data.frame(
   sample=training_samples,
   cancer_type=getMetadata(training_samples, 'cancer_type'),
   cohort=getMetadata(training_samples, 'cohort')
)

ct_counts <- table(cancer_type=training_samples$cancer_type, cohort=training_samples$cohort)
ct_counts <- as.data.frame(unclass(ct_counts))
training_samples$is_selected.HMF <-
   training_samples$cohort=='HMF' &
   training_samples$cancer_type %in% rownames(ct_counts)[ct_counts$HMF>=15]

training_samples$is_selected.PCAWG <-
   training_samples$cohort=='PCAWG' &
   training_samples$cancer_type %in% rownames(ct_counts)[ct_counts$PCAWG>=15]

write.table(
   training_samples,
   paste0(out_dir,'/training_sample_selection.txt'),
   sep='\t',quote=F, row.names=F
)

## Write features ================================
## HMF
message('Writing features; HMF training samples...')
(function(){
   sel_samples <- training_samples$sample[training_samples$is_selected.HMF]
   df <- df_features[sel_samples,]

   write.table(
      df,
      gzfile(paste0(out_dir,'/features.HMF.txt.gz')),
      sep='\t',quote=F
   )

   write.table(
      data.frame(sample=rownames(df), response=df$response),
      paste0(out_dir,'/training_labels.HMF.txt'),
      sep='\t',quote=F
   )
})()

message('Writing features; PCAWG training samples...')
(function(){
   sel_samples <- training_samples$sample[training_samples$is_selected.PCAWG]
   df <- df_features[sel_samples,]

   write.table(
      df,
      gzfile(paste0(out_dir,'/features.PCAWG.txt.gz')),
      sep='\t',quote=F
   )

   write.table(
      data.frame(sample=rownames(df), response=df$response),
      paste0(out_dir,'/training_labels.PCAWG.txt'),
      sep='\t',quote=F
   )
})()

message('Writing features; all samples...')
(function(){
   write.table(
      df_features,
      gzfile(paste0(out_dir,'/features_allSamples.txt.gz')),
      sep='\t',quote=F
   )
})()


