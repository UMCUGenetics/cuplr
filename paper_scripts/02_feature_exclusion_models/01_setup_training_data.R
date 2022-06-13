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

model_dir <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21c_featGroupExcl/')
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

## Select training samples
sample_whitelist <- metadata[metadata$is_training_sample,'sample_id']
training_samples <- rownames(df_features)[rownames(df_features) %in% sample_whitelist]

#View(metadata[metadata$in_train_set,])

## Export ================================
## Training labels
message('Writing labels...')
training_labels <- data.frame(
   sample=rownames(df_features),
   response=df_features$response,
   in_training_set=rownames(df_features) %in% training_samples
)
training_labels <- rbind(
   training_labels[match(training_samples, training_labels$sample),],
   subset(training_labels, !in_training_set)
)

write.table(
   training_labels,
   paste0(out_dir,'/training_labels.txt'),
   sep='\t',quote=F,row.names=F
)

## Feature names
write.table(
   matrix(colnames(df_features[,-1]), ncol=1),
   paste0(out_dir,'/feature_names.txt'),
   sep='\t',quote=F,row.names=F,col.names=F
)

## Features
message('Writing features; training samples...')
(function(){

   ## Subset features
   df <- df_features[training_samples,]
   df$response <- droplevels(df$response)

   write.table(
      df,
      gzfile(paste0(out_dir,'/features.txt.gz')),
      sep='\t',quote=F
   )
   #saveRDS(df, paste0(out_dir,'/features.rds'))
})()

message('Writing features; all samples...')
(function(){
   write.table(
      df_features,
      gzfile(paste0(out_dir,'/features_allSamples.txt.gz')),
      sep='\t',quote=F
   )
   #saveRDS(df_features, paste0(out_dir,'/features_allSamples.rds'))
})()
