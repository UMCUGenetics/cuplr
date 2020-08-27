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
#' Perform cross validation
#'
#' @param df A dataframe containing the training features and response column.
#' @param folds Alternative input to df. A list with the structure: `l[[1]][c('train','test')]`
#' where 'train' and 'test' are dataframes which features and a response column. This is the output
#' from `createCvTrainTestSets()`
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
   df=NULL, folds=NULL, train.func, train.func.args=list(), colname.response='response',
   k=10, seed=1, multi.core=F, verbose=1
){
   #df=training_data
   #colname.response='response'
   #k=5
   #verbose=2
   #train.func=trainNeuralNet
   #train.func.args=list(feature.groups=feature_groups)

   if(!is.list(train.func.args)){ stop('`train.func.args` must be a list') }

   if(!is.null(df)){
      folds <- mltoolkit::createCvTrainTestSets(df, stratify.by.col=colname.response, k=k)
   }
   if(is.null(folds)){ stop('Either `df` or `folds` must be provided') }

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

####################################################################################################
#' Train a random forest (for binary classification)
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
#' @param n.resamples.true Number of resampling values to generate for the TRUE class
#' @param n.resamples.false Number of resampling values to generate for the FALSE class
#' @param min.size.diff Default=30. Minimum difference between resampling values (integer).
#' @param midpoint.type Can be 'geometric', 'arithmetic', or 'none'. Calculate the resampling values
#' from a->midpoint and b->midpoint? If 'none', resampling values will be calculated from a->b and
#' b->a.
#' @param max.upsample.ratio Default=10. Remove pairs where a or b are upsampled more than this
#' value
#' @param ntree Number of decision trees in the random forest
#' @param get.local.increments If TRUE, get local increments to calculate feature contributions
#' @param calc.imp Calculate feature importance?
#' @param seed Random seed as an integer
#' @param fold.num Cross validation fold number. Only used in the prefix of progress messages. If
#' not specified, no prefix will be displayed.
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForest <- function(
   train, test=NULL, colname.response='response',

   ## univariate feature slection
   do.feat.sel=TRUE, feat.sel.max.qvalue=0.01, feat.sel.max.pvalue=1, feat.sel.top.n.features=NULL,

   ## class balancing
   balance.classes=F, k.inner=5,
   n.resamples.true=4, n.resamples.false=4, midpoint.type='geometric',
   min.size.diff=NULL, max.upsample.ratio=10,

   ntree=500, get.local.increments=TRUE,
   calc.imp=T,
   seed=NULL, fold.num=NULL, verbose=1
){
   if(F){
      train=folds[[1]]$train
      test=folds[[1]]$test
      colname.response='response'

      ## univariate feature slection
      do.feat.sel=TRUE; feat.sel.max.qvalue=0.01; feat.sel.max.pvalue=1; feat.sel.top.n.features=NULL;

      ## class balancing
      balance.classes=F; k.inner=5;
      n.resamples.true=4; n.resamples.false=4; midpoint.type='geometric';
      min.size.diff=NULL; max.upsample.ratio=10

      ntree=500; get.local.increments=TRUE;
      calc.imp=T;
      seed=NULL; fold.num=NULL; verbose=1
   }

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
   y <- train_data$y

   if(do.feat.sel){
      if(verbose){ message(msg_prefix,'> Performing feature selection...') }
      x <- univarFeatSel(
         train_data$x, y,
         max.qvalue=feat.sel.max.qvalue, max.pvalue=feat.sel.max.pvalue,
         sel.top.n.features=feat.sel.top.n.features,
         verbose=(verbose==3)
      )
   } else {
      x <- train_data$x
   }

   ##----------------------------------------------------------------------
   out <- list()
   out$model <- NA ## placeholder

   if(balance.classes){
      if(verbose){ message(msg_prefix,'> Balancing classes...') }
      resampling_grid <- resamplingGrid(
         a=sum(y==TRUE), b=sum(y==FALSE),
         breaks.a=n.resamples.true, breaks.b=n.resamples.false,
         min.size.diff=min.size.diff, midpoint.type=midpoint.type, max.upsample.ratio=max.upsample.ratio
      )

      train_new <- structure(cbind(y, x), names=c(colname.response, colnames(x)))
      perfs_inner_cv <- lapply(1:nrow(resampling_grid), function(i){
         #i=6
         size_true <- resampling_grid[i,'size_a']
         size_false <- resampling_grid[i,'size_b']

         if(verbose==2){
            ratio_true <- round(resampling_grid[i,'ratio_a'], 2)
            ratio_false <- round(resampling_grid[i,'ratio_b'], 2)
            message(msg_prefix,'>> TRUE, FALSE: ',ratio_true,'x, ', ratio_false,'x')
         }

         inner_folds <- mltoolkit::createCvTrainTestSets(train_new, stratify.by.col=colname.response, k=k.inner)

         for(j in 1:length(inner_folds)){
            inner_folds[[j]]$train <- resampleClasses(
               df=inner_folds[[j]]$train, colname.response=colname.response,
               target.sample.sizes=c('TRUE'=size_true,'FALSE'=size_false)
            )
         }

         inner_cv <- crossValidate(
            folds=inner_folds, colname.response=colname.response,
            train.func=trainRandomForest,
            train.func.args=list(ntree=ntree, do.feat.sel=F, get.local.increments=F, calc.imp=F),
            verbose=(verbose==3)
         )

         unlist(lapply(inner_cv, function(j){
            #j=inner_cv[[1]]
            pred <- cbind(
               as.data.frame(j$test_set$probabilities),
               response=j$test_set$actual
            )

            confusion <- mltoolkit::confusionMatrix(
               predicted=pred[,'TRUE'],
               actual=pred[,'response']
            )

            perf <- mltoolkit::calcPerfCompound(confusion, compound.metric='pr', metric.names.as.x.y=T)
            mltoolkit::calcAUC(x=perf[,'x'], y=perf[,'y'])
         }))
      })

      resampling_grid$perf <- sapply(perfs_inner_cv, mean)
      resampling_grid$sd <- sapply(perfs_inner_cv, sd)
      resampling_grid$best <- resampling_grid$perf==max(resampling_grid$perf)

      # resampling_grid$index <- 1:nrow(resampling_grid)
      # ggplot(resampling_grid, aes(x=index, y=perf)) +
      #    geom_bar(stat='identity') +
      #    geom_linerange(aes(ymin=perf-sd, ymax=perf+sd))

      ## Apply class resampling
      train_new <- resampleClasses(
         df=train_new, colname.response=colname.response,
         target.sample.sizes=
            structure(
               unlist(subset(resampling_grid, best, c(size_a,size_b))),
               names=c('TRUE','FALSE')
            )
      )
      x <- train_new[,colnames(train_new)!=colname.response]
      y <- train_new[,colname.response]
      rm(train_new)

      out$resampling_grid <- resampling_grid
   }

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'> Training random forest...') }
   model <- randomForest::randomForest(
      x=x,
      y=if(is.logical(y)){ factor(y, c('TRUE','FALSE')) } else { y }, strata=y,
      proximity=F, ntree=ntree, importance=T,
      keep.inbag=T, replace=F, ## required for calculating local increments
      na.action=na.roughfix,
      do.trace=(verbose==3)
   )

   if(get.local.increments){
      if(verbose){ message(msg_prefix,'> Getting local increments...') }
      ## Used for calculating per sample feature contributions
      invisible(capture.output(
         model$localIncrements <- rfFC::getLocalIncrements(model, x)
      ))
   }

   out$model <- model
   out$categorical_lvls <- categorical_lvls

   ##----------------------------------------------------------------------
   if(calc.imp){
      if(verbose){ message(msg_prefix,'> Calculating feature importance...') }
      out$imp <- t(randomForest::importance(model, type=1))
   }

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){
         message(msg_prefix,'> Predicting on test set...')
      }
      out$test_set <- list(
         probabilities = randomForest:::predict.randomForest(model, test_data$x, type='prob'),
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

   out$seed <- seed

   return(out)
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
   calc.imp=T, ntree=500, get.local.increments=TRUE,
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

      out$imp <- do.call(rbind, lapply(model, function(i){
         #i=model[[1]]
         df <- randomForest::importance(i, type=1)
         v <- structure(df[,1],names=rownames(df))
         v <- v[exist_features]
         names(v) <- exist_features
         v[is.na(v)] <- 0
         return(v)
      }))
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





