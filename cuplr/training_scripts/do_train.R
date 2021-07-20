options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')


library(randomForest)
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/statsExtra/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/nmf/'))
#devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))

## Main ================================
args <- commandArgs(trailingOnly=T)

if(F){
   wd <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.16d_newBlacklist/')
   training_data_path <- paste0(wd,'/features/features.rds')
   fold_indexes <- paste0(wd,'/cv_out/01/fold_indexes.rds')
   out_path <- paste0(wd,'/cv_out/01/model.rds')
   seed <- 1
}

## --------------------------------
message('## Loading training data...')
training_data_path <- args[1] ## <<<<<<<<
print(training_data_path)

training_data <- readRDS(training_data_path)
training_data$response <- as.factor(training_data$response)

#is_character <- sapply(training_data, is.character)
#colnames(training_data)[is_character]
#head2(training_data[,is_character])
#v <- training_data[,'rmd.6p_47']
#v[!grepl('(0|[.])',v)]

## --------------------------------
message('## Loading folding indexes...')
fold_indexes <- args[2] ## <<<<<<<<
print(fold_indexes)
fold_indexes <- if(is.na(fold_indexes) || nchar(fold_indexes)==0){
   list()
} else {
   readRDS(fold_indexes)
}

## --------------------------------
out_path <- args[3] ## <<<<<<<<

seed <- args[4] ## <<<<<<<<
if(is.na(seed)){ seed <- 1 }

## Data --------------------------------
message('## Subsetting data...')
if(length(fold_indexes)!=0){
   train <- training_data[fold_indexes$train,]
   test <- training_data[fold_indexes$test,]
} else {
   train <- training_data
   test <- NULL
}
rm(training_data)

## Train --------------------------------
feature_names <- colnames(train)[colnames(train)!='response']
alternative <- structure(
   rep('greater',length(feature_names)),
   names=feature_names
)
alternative[grep('^rmd',names(alternative))] <- 'two.sided'
alternative[grep('^mut_load',names(alternative))] <- 'two.sided'
alternative[names(alternative)=='sv.n_events'] <- 'two.sided'
alternative[names(alternative)=='gender.gender'] <- 'two.sided'

train_func_args <- list(
   #inner.holdout.fraction=c(1,3),
   do.rmd.nmf=T,
   #tmp.dir=paste0(dirname(out_path),'/tmp/')
   args.trainRandomForest=list(
      ntree=500,
      do.feat.sel=T,
      feat.sel.alternative=alternative,
      feat.sel.whitelist='gender.gender',
      feat.sel.max.pvalue=0.01,
      feat.sel.min.cliff.delta=0.1,
      feat.sel.min.cramer.v=0.1,
      feat.sel.top.n.features=100,
      balance.classes='resample',
      verbose=2
   ),
   tmp.dir=paste0(dirname(out_path),'/tmp/'),
   rm.tmp.dir=F
)

if(!file.exists(out_path)){
   message('\n## Training model...')
   model <- do.call(
      trainRandomForestEnsemble,
      c(
         list(train=train, seed=seed, multi.core=T, verbose=2),
         train_func_args
      )
   )
   saveRDS(model, out_path)
} else {
   message('\n## Loading model...')
   model <- readRDS(out_path)
}

if(!is.null(test)){
   message('\n## Making test set report...')
   test_set_report <- predict.randomForestEnsemble(
      object=model, newdata=test,
      type='report', verbose=T
   )

   test_set_report$class_actual <- test$response
   test_set_report$imp <- model$imp

   saveRDS(
      test_set_report,
      paste0(dirname(out_path),'/test_set_report.rds')
   )

}

