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
      paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/models/0.06c_probWeighRf_balanceClasses/features/features.rds')
   )
   training_data$response <- training_data$response=='Lymphoid'

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

   set.seed(seed)
   if(!is.null(df)){
      folds <- mltoolkit::createCvTrainTestSets(df, stratify.by.col=colname.response, k=k)
   }
   if(is.null(folds)){ stop('Either `df` or `folds` must be provided') }

   if(verbose){ message('Performing CV...') }
   main <- function(i){
      #i=1
      #if(verbose){ message('\n## Fold: ',i) }
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
   do.feat.sel=T, feat.sel.max.qvalue=NULL, feat.sel.max.pvalue=0.01, feat.sel.top.n.features=NULL,

   ## class balancing
   balance.classes=F, k.inner=5,
   n.resamples.true=4, n.resamples.false=4, midpoint.type='geometric',
   min.size.diff=NULL, max.upsample.ratio=10,

   ntree=500, get.local.increments=T,
   seed=NULL, fold.num=NULL, verbose=1
){
   if(F){
      train=folds[[1]]$train
      test=folds[[1]]$test
      colname.response='response'

      ## univariate feature slection
      do.feat.sel=T; feat.sel.max.qvalue=NULL; feat.sel.max.pvalue=0.01; feat.sel.top.n.features=25;

      ## class balancing
      balance.classes=F; k.inner=5;
      n.resamples.true=4; n.resamples.false=4; midpoint.type='geometric';
      min.size.diff=NULL; max.upsample.ratio=10

      ntree=500; get.local.increments=T;
      calc.imp=T;
      seed=NULL; fold.num=NULL; verbose=2
   }

   if(!is.null(seed)){ set.seed(seed) }

   msg_prefix <- if(!is.null(fold.num)){ paste0('[[',fold.num,']] ') } else { '' }

   ##----------------------------------------------------------------------
   if(verbose>=1){message(msg_prefix,'> Preparing input features and response variable...')}
   train_data <- dfToFeaturesAndResponse(train, colname.response=colname.response)
   
   if(any(sapply(train_data$x, is.character))){
      stop('Categorical features must be factors and not characters')
   }
   
   if(!is.null(test)){
      test_data <- dfToFeaturesAndResponse(test, colname.response=colname.response)
      if(any(sapply(test_data$x, is.character))){
         stop('Categorical features must be factors and not characters')
      }
   }
   
   ##----------------------------------------------------------------------
   y <- train_data$y
   #y <- y=='Stomach'

   x <- train_data$x

   if(do.feat.sel){
      if(verbose>=1){ message(msg_prefix,'> Performing feature selection...') }
      x <- univarFeatSel(
         x, y,
         max.qvalue=feat.sel.max.qvalue, max.pvalue=feat.sel.max.pvalue,
         sel.top.n.features=feat.sel.top.n.features,
         verbose=(verbose>=3)
      )
   }

   ##----------------------------------------------------------------------
   out <- list()
   out$model <- NA ## placeholder

   if(balance.classes){
      if(verbose>=1){ message(msg_prefix,'> Balancing classes...') }
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

         if(verbose>=2){
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
            train.func.args=list(ntree=ntree, do.feat.sel=F, get.local.increments=F),
            verbose=(verbose>=3)
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
      resampling_grid$best <- FALSE
      resampling_grid$best[ which.max(resampling_grid$perf) ] <- TRUE

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
   if(verbose>=1){ message(msg_prefix,'> Training random forest...') }
   model <- randomForest::randomForest(
      x=x,
      y=if(is.logical(y)){ factor(y, c('TRUE','FALSE')) } else { y }, strata=y,
      proximity=F, ntree=ntree, importance=T,
      keep.inbag=T, replace=F, ## required for calculating local increments
      na.action=na.roughfix,
      do.trace=F
   )

   if(get.local.increments){
      if(verbose>=1){ message(msg_prefix,'> Getting local increments...') }
      ## Used for calculating per sample feature contributions
      invisible(capture.output(
         model$localIncrements <- rfFC::getLocalIncrements(model, x)
      ))
   }

   out$model <- model

   ## Store levels from categorical features
   categorical_lvls <- lapply(x, levels)
   categorical_lvls <- categorical_lvls[ !sapply(categorical_lvls, is.null) ]
   categorical_lvls <- categorical_lvls[names(categorical_lvls) %in% colnames(x)]

   out$categorical_lvls <- categorical_lvls

   ##----------------------------------------------------------------------
   if(verbose>=1){ message(msg_prefix,'> Calculating feature importance...') }
   out$imp <- sort(
      randomForest::importance(model, type=1)[,1],
      decreasing=T
   )

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose>=1){
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
      # df <- cbind(
      #    as.data.frame(out$test_set$probabilities),
      #    class=out$test_set$actual,
      #    response=out$test_set$actual=='Stomach'
      # )
      # df[order(df$response, decreasing=T),]
      # df[order(df[,'TRUE'], decreasing=T),]
   }

   out$seed <- seed
   class(out) <- 'randomForestContainer'

   return(out)
}

print.randomForestContainer <- function(rf){
   cat(
      'Objects in list:\n',
      paste0('$',names(rf))
   )

   cat('\n\n$model'); print(rf$model)

   cat(
      '\n$categorical_lvls',
      '\nNumber of categorical features:',
      length(rf$categorical_lvls),
      '\n'
   )

   cat(
      '\n$imp',
      '\nNumber of features:', length(rf$imp),
      '\nTop features:\n'
   )
   print(head(rf$imp, 5))
}

####################################################################################################
#' Train a random forest ensemble (for multiclass classification)
#'
#' @param train A dataframe containing the training features and response column.
#' @param test A dataframe containing the training features and response column. Used for testing
#' the performance of the trained model
#' @param colname.response The column name of the response variable (i.e. training labels)
#' @param args.trainRandomForest A list of arguments that can be passed to `trainRandomForest()`
#' @param seed Random seed as an integer
#' @param fold.num Cross validation fold number. Only used in the prefix of progress messages. If
#' not specified, no prefix will be displayed.
#' @param multi.core Use multiple cores?
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForestEnsemble <- function(
   train, test=NULL, colname.response='response',
   args.trainRandomForest=list(),
   seed=NULL, fold.num=NULL, multi.core=F, verbose=1,
   tmp.dir=NULL
){
   if(F){
      train=folds[[1]]$train
      test=folds[[1]]$test
      colname.response='response'

      args.trainRandomForest=list(feat.sel.max.pvalue=0.01, feat.sel.top.n.features=100)

      seed=NULL
      fold.num=1
      verbose=2
   }

   if(!is.null(seed)){ set.seed(seed) }

   msg_prefix <- if(!is.null(fold.num)){ paste0('[[',fold.num,']] ') } else { '' }

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'> Training random forest ensemble...') }

   y_ohe <- oneHotEncode(as.factor(train[,colname.response]))
   
   main <- function(i){
      #i=16
      #i='Lymphoid'
      
      if(verbose==2){ message(msg_prefix, '[',i,'/',ncol(y_ohe),']: ', colnames(y_ohe)[i] ) }
      
      y <- unname(y_ohe[,i])
      train_new <- train
      train_new$response <- y_ohe[,i]
      
      rf <- do.call(
         trainRandomForest,
         c(list(train=train_new), args.trainRandomForest)
      )
      
      return(rf)
   }
   
   if(!multi.core){
      ensemble <- lapply(1:ncol(y_ohe), main)
   } else {
      require(foreach)
      require(parallel)
      require(doParallel)
      
      #cores <- Sys.getenv("SLURM_NTASKS_PER_NODE")
      cores <- parallel::detectCores()
      if(verbose){ message('Multicore: ',cores,' cores') }
      
      registerDoParallel(cores=cores)
      ensemble <- foreach(i=1:ncol(y_ohe)) %dopar% { main(i) }
      #stopCluster(cl)
   }
   names(ensemble) <- colnames(y_ohe)
   
   if(!is.null(tmp.dir)){
      dir.create(tmp.dir, showWarnings=F)
      saveRDS(ensemble, paste0(tmp.dir,'/model.rds'))
   }
   #ensemble=readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.06e_balanceClasses/tmp/model.rds')
   
   out <- list()
   class(out) <- c('list','randomForestEnsemble')
   
   out$ensemble <- ensemble
   out$seed <- seed

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'> Calculating feature importance...') }
   exist_features <- unique(unlist(lapply(ensemble, function(i){
      names(i$model$forest$ncat)
   })))
   features <- colnames(train)[!(colnames(train) %in% colname.response)]
   exist_features <- features[ features %in% exist_features ] ## Preserve original feature order
   rm(features)

   out$imp <- do.call(rbind, lapply(ensemble, function(i){
      #i=ensemble[[1]]
      df <- randomForest::importance(i$model, type=1)
      v <- structure(df[,1],names=rownames(df))
      v <- v[exist_features]
      names(v) <- exist_features
      v[is.na(v)] <- 0
      return(v)
   }))

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){
         message(msg_prefix,'> Predicting on test set...')
      }
      out$test_set <- list(
         probabilities = predict.randomForestEnsemble(
            out, 
            test[,!(colnames(test) %in% colname.response)], 
            type='prob'
         ),
         predicted = NA,
         #actual = test[,colname.response]
         actual = factor(test[,colname.response], colnames(y_ohe))
      )

      out$test_set$predicted <- with(out$test_set,{
         factor(
            colnames(probabilities)[max.col(probabilities)],
            levels=colnames(y_ohe)
         )
      })
   }

   return(out)
}

print.randomForestEnsemble <- function(object){
   cat('$ensemble\n\nBinary RF names:\n')
   print(names(object$ensemble))
   
   cat('\nOther list levels:\n')
   cat( paste0('$',names(object)[3:length(object)]) )
}

#' Predict method for random forest ensemble
#'
#' @param object An object of class randomForestEnsemble
#' @param newdata A data frame or matrix containing new data. (Note: If not given, the out-of-bag
#' prediction in object is returned)
#' @param type 'response', 'prob' or 'votes', indicating the type of output: predicted values,
#' matrix of class probabilities, or matrix of vote counts. 'class' is allowed, but automatically
#' converted to "response", for backward compatibility.
#'
#' @return 'response': predicted classes (the classes with majority vote). 'prob' matrix of class
#' probabilities (one column for each class and one row for each input). 'vote'
#' matrix of vote counts (one column for each class and one row for each new input); either in raw
#' counts or in fractions (if norm.votes=TRUE).
#' @export
#'
predict.randomForestEnsemble <- function(object, newdata, type='response', verbose=F){
   # object=ensemble
   # newdata=test[,!(colnames(test) %in% colname.response)]
   # type='prob'
   
   ensemble <- object$ensemble
   
   ## Force factor levels from training set onto new data
   categorical_feature_names <- names(ensemble[[1]]$categorical_lvls)
   x <- as.data.frame(lapply(colnames(newdata), function(i){
      #i='gene_def.AR'
      if(!(i %in% categorical_feature_names)){
         return(newdata[,i])
      } else{
         factor(newdata[,i], levels=ensemble[[1]]$categorical_lvls[[i]])
      }
   }))
   rownames(x) <- rownames(newdata)
   colnames(x) <- colnames(newdata)
   
   if(verbose){ counter <- 0 }
   m <- do.call(cbind, lapply(ensemble, function(i){
      #i=ensemble$Breast
      if(verbose){
         counter <<- counter + 1
         message('[RF ',counter,'/',length(ensemble),']: ', names(ensemble)[counter])
      }
      randomForest:::predict.randomForest(i$model, x, type='prob')[,1]
   }))
   
   if(type=='prob'){
      m
   } else {
      factor( colnames(m)[ max.col(m) ], levels=colnames(m) )
   }
}




