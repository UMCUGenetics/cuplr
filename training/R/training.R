## Data for debugging ================================
if(F){
   base_dir <- list(
      hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
      mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
      local='/Users/lnguyen/Documents/projects/P0013_WGS_patterns_Diagn/'
   )

   for(i in base_dir){
      if(dir.exists(i)){
         base_dir <- i
         break
      }
   }

   devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/'))
   library(mltoolkit)

   # training_data <- read.delim(
   #    paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/models/0.06a_ownDrivers/features.txt.gz'),
   #    check.names=F, stringsAsFactors=T
   # )

   training_data <- readRDS(
      paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/models/0.06a_ownDrivers/features/features.rds')
   )

   folds <- createCvTrainTestSets(training_data)

}


####################################################################################################
#' Train a random forest ensemble (for multiclass classification)
#'
#' @param train A dataframe containing the training features and response column.
#' @param test A dataframe containing the training features and response column. Used for testing
#' the performance of the trained model
#' @param colname.response The column name of the response variable (i.e. training labels)
#' @param do.feat.sel Perform univariate feature selection?
#' @param feat.sel.max.qvalue Only features with q-value (from univariate feature selection) lower
#' than this will be kept
#' @param feat.sel.max.pvalue Only features with p-value (from univariate feature selection) lower
#' than this will be kept
#' @param feat.sel.top.n.features Only keep the top number of features
#' @param calc.imp Calculate feature importance?
#' @param imp.metric Which type of feature importance to calculate
#' @param ntree Number of decision trees in the random forest
#' @param get.local.increments If TRUE, get local increments to calculate feature contributions
#' @param seed Random seed as an integer
#' @param fold.num Cross validation fold number. Only used in the prefix of progress messages. If
#' not specified, no prefix will be displayed.
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForestEnsemble <- function(
   train, test=NULL, colname.response='response',
   do.feat.sel=TRUE, feat.sel.max.qvalue=0.01, feat.sel.max.pvalue=1, feat.sel.top.n.features=NULL,
   calc.imp=T, imp.metric='mda',
   ntree=500, get.local.increments=TRUE,
   seed=NULL, fold.num=NULL, verbose=1
){
   # train=folds[[1]]$train
   # test=folds[[1]]$test
   # colname.response='response'
   # do.feat.sel=TRUE
   # feat.sel.max.qvalue=1
   # feat.sel.max.pvalue=0.01
   # feat.sel.top.n.features=200
   # calc.imp=T
   # imp.metric='mda'
   # ntree=200
   # get.local.increments=T
   # seed=NULL
   # fold.num=1
   # verbose=2

   if(!is.null(seed)){ set.seed(seed) }

   msg_prefix <- if(!is.null(fold.num)){ paste0('[[',fold.num,']] ') } else { '' }

   ##----------------------------------------------------------------------
   if(verbose){message(msg_prefix,'> Preparing input features and response variable...')}
   train_data <- dfToFeaturesAndResponse(train, colname.response=colname.response)

   if(any(sapply(train_data$x, is.character))){
      stop('Categorical features must be factors and not characters')
   }
   categorical_lvls <- lapply(train_data$x, levels)
   categorical_lvls <- categorical_lvls[ !sapply(categorical_lvls, is.null) ]

   if(!is.null(test)){
      test_data <- dfToFeaturesAndResponse(test, colname.response=colname.response)
   }

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'> Training random forest ensemble...') }
   model <- lapply(1:ncol(train_data$y_ohe),function(i){
      #i='Prostate'

      if(verbose==2){ message(msg_prefix, '[',i,'/',ncol(train_data$y_ohe),']: ', colnames(train_data$y_ohe)[i] ) }

      y <- unname(train_data$y_ohe[,i])

      if(do.feat.sel){
         if(verbose==2){ message(msg_prefix,'>> Performing feature selection...') }
         x <- univarFeatSel(
            train_data$x, y,
            max.qvalue=feat.sel.max.qvalue, max.pvalue=feat.sel.max.pvalue,
            sel.top.n.features=feat.sel.top.n.features,
            verbose=F
         )
      } else {
         x <- train_data$x
      }

      if(verbose==2){ message(msg_prefix,'>> Training random forest...') }
      rf <- randomForest::randomForest(
         x=x,
         y=factor(y, levels=c('TRUE','FALSE')),
         strata=y, proximity=F, ntree=ntree,
         importance=(calc.imp & imp.metric=='mda'),
         keep.inbag=T, replace=F, ## required for calculating local increments
         na.action=na.roughfix,
         do.trace=F
      )

      if(get.local.increments){
         if(verbose==2){ message(msg_prefix,'>> Getting local increments...') }
         ## Used for calculating per sample feature contributions
         invisible(capture.output(
            rf$localIncrements <- rfFC::getLocalIncrements(rf, x)
         ))
      }

      return(rf)
   })
   names(model) <- colnames(train_data$y_ohe)
   class(model) <- c('list','randomForestEnsemble')

   out <- list()
   out$model <- model
   out$seed <- seed

   ##----------------------------------------------------------------------
   if(calc.imp){
      if(verbose){ message(msg_prefix,'> Calculating feature importance...') }
      exist_features <- unique(unlist(lapply(model, function(i){
         names(i$forest$ncat)
      })))
      exist_features <- colnames(train_data$x)[ colnames(train_data$x) %in% exist_features ] ## Preserve original feature order

      if(imp.metric=='mda'){
         out$imp <- do.call(rbind, lapply(model, function(i){
            #i=model[[1]]
            df <- randomForest::importance(i, type=1)
            v <- structure(df[,1],names=rownames(df))
            v <- v[exist_features]
            names(v) <- exist_features
            v[is.na(v)] <- 0
            return(v)
         }))

      } else {
         if(!is.null(test)){
            out$imp <- permutationImportance(
               model, test_data$x[,exist_features], test_data$y, metric=imp.metric,
               verbose=(verbose==2)
            )
         } else {
            stop(msg_prefix,'> Feature importance calculation requires test data when using the permutation method')
         }
      }
   }

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){
         message(msg_prefix,'> Predicting on test set...')
      }
      out$test_set <- list(
         probabilities = predict.randomForestEnsemble(model, test_data$x, type='prob'),
         predicted = NA,
         actual = test_data$y
      )

      out$test_set$predicted <- with(out$test_set,{
         factor(
            colnames(probabilities)[max.col(probabilities)],
            levels=colnames(probabilities)
         )
      })
   }

   return(out)
}

####################################################################################################
#' Perform cross validation
#'
#' @param df A dataframe containing the training features and response column.
#' @param train.func A training function.
#' @param train.func.args Arguments (as a list) that can be passed to the training function.
#' @param colname.response The column name of the response variable (i.e. training labels).
#' @param k Number of cross validation folds.
#' @param seed Random seed as in integer.
#' @param multi.core Use multiple cores?
#' @param verbose Show progress messages?
#'
#' @return A list containing for each CV fold, the output of the training function
#' @export
#'
crossValidate <- function(
   df, train.func, train.func.args=list(), colname.response='response',
   k=10, seed=1, multi.core=F, verbose=1
){
   #df=training_data
   #colname.response='response'
   #k=5
   #verbose=2
   #train.func=trainNeuralNet
   #train.func.args=list(feature.groups=feature_groups)

   if(!is.list(train.func.args)){ stop('`train.func.args` must be a list') }

   folds <- createCvTrainTestSets(df, stratify.by.col=colname.response, k=k)

   if(verbose){ message('Performing CV...') }
   main <- function(i){
      #i=1
      if(verbose){ message('\n## Fold: ',i) }
      fold <- folds[[i]]

      fold_seed <- as.integer(paste0(seed,i))
      set.seed(fold_seed)

      out <- do.call(
         train.func,
         c(
            list(train=fold$train, test=fold$test, fold.num=i, verbose=verbose),
            train.func.args
         )
      )

      out$seed <- fold_seed
      return(out)
   }

   #out <- main(1)
   if(!multi.core){
      lapply(1:length(folds), main)
   } else {
      require(foreach)
      require(parallel)
      require(doParallel)

      #cores <- Sys.getenv("SLURM_NTASKS_PER_NODE")
      cores <- parallel::detectCores()
      if(verbose){ message('Multicore: ',cores,' cores') }

      registerDoParallel(cores=cores)
      foreach(i=1:length(folds)) %dopar% { main(i) }
      #stopCluster(cl)
   }
}

# cv_out <- crossValidate(
#    training_data, train.func=trainNeuralNet,
#    train.func.args=list(feature.groups=feature_groups, calc.imp=T),
#    k=3, verbose=2
# )
#
# cv_out <- crossValidate(
#    training_data, train.func=trainRandomForest,
#    train.func.args=list(calc.imp=T),
#    k=3, verbose=2
# )
# cv_out[[1]]$imp
