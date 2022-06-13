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
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))

## Main ================================================================
wd <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21b_DR104-update5_pcawgSoftFilt2_clonal/')

training_labels <- read.delim(paste0(wd,'/features/training_labels.txt'))
training_labels <- subset(training_labels, in_training_set)

spawnCvJobs(
   df=training_labels,
   train.data.path=paste0(wd,'/features/features.txt.gz'),
   train.script.path=paste0(wd,'/do_train.R'),
   k=15, time='4:00:00', mem='60G', n.tasks.per.node=20, seed=5
)
