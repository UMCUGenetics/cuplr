options(stringsAsFactors=F)

## Paths ================================================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')

devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))

## Main ================================================================
wd <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21c_featGroupExcl/')

##
training_labels <- read.delim(paste0(wd,'/features/training_labels.txt'))
training_labels <- subset(training_labels, in_training_set)

##
feature_names <- read.delim(paste0(wd,'/features/feature_names.txt'), header=F)
feature_names <- feature_names[,1]
feature_groups <- groupFeaturesByTag(feature_names)

##
seeds <- structure(rep(17, length(feature_groups)), names=names(feature_groups))
cv_out_parent_dir <- paste0(wd,'/seed_models/17/cv_out/')
dir.create(cv_out_parent_dir, showWarnings=F, recursive=T)

counter <- 0
for(i in names(feature_groups)){
   #i=names(feature_groups)[1]
   counter <- counter + 1

   counter_f <- formatC(counter, width=2, format="d", flag="0")
   cv_out_dir <- paste0(cv_out_parent_dir,'/',counter_f,'_',i,'/')

   if(dir.exists(cv_out_dir)){
      message('Skipping: ',i,'. Dir exists: ',cv_out_dir)
      next
   }

   dir.create(cv_out_dir, showWarnings=F)

   ##
   feat_excl_path <- paste0(cv_out_dir,'/feat_excl.rds')
   saveRDS(
      feature_groups[[i]],
      feat_excl_path
   )

   seed <- seeds[i]

   spawnCvJobs(
      df=training_labels,
      train.data.path=paste0(wd,'/features/features.txt.gz'),
      train.script.path=paste0(wd,'/do_train.R'),
      trailing.args=feat_excl_path,
      cv.out.dir=cv_out_dir,
      k=15, time='4:00:00', mem='60G', n.tasks.per.node=20, seed=seed
   )

  #if(counter>=4){ break }
}


