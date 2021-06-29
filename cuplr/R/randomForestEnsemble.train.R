# ## Data for debugging ================================
# if(F){
#    base_dir <- list(
#       hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
#       mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
#       local='/Users/lnguyen/Documents/projects/P0013_WGS_patterns_Diagn/'
#    )
#
#    for(i in base_dir){
#       if(dir.exists(i)){
#          base_dir <- i
#          break
#       }
#    }
#
#    devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/commonUtils/'))
#    devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
#    library(mltoolkit)
#
#    features <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10f_rmd1mbSigs/features/features.txt.gz')
#    features$response <- as.factor(features$response)
#    fold_indexes <- readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10f_rmd1mbSigs/cv_out/01/fold_indexes.rds')
# }

####################################################################################################
#' Remove columns by name
#'
#' @param df A dataframe or matrix
#' @param columns Column names
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#'
#' @return A dataframe or matrix
#'
rmColumns <- function(df, columns, drop=F){
   if(!is.character(columns)){ stop('`columns must be a character vector`') }
   df[,!(colnames(df) %in% columns), drop=drop]
}

#' Prepares a feature dataframe for multiclass classifier training
#'
#' @description Separates a dataframe containing features and response variable. The response
#' variable is also one hot encoded
#'
#' @param df A dataframe
#' @param colname.response Name of the response column
#'
#' @return A list containing the feature matrix/dataframe, response vector, and response
#' one-hot encode matrix
#' @export
#'
dfToFeaturesAndResponse <- function(df, colname.response='response'){
   x <- df[,colnames(df)!=colname.response,drop=F]
   #x <- as.matrix(x)

   if(is.logical(df[,colname.response])){
      y <- df[,colname.response]
      y_ohe <- NA
   } else {
      y <- as.factor(df[,colname.response])
      y_ohe <- oneHotEncode(y, sample.names=rownames(x))
   }

   list(x=x, y=y, y_ohe=y_ohe)
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

      cv_msg_prefixt <- paste0('[Fold: ',i,']')
      out <- do.call(
         train.func,
         c(
            list(train=fold$train, test=fold$test, msg.prefix=i, verbose=verbose),
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
#' @param feat.sel.alternative Can be 'two.sided', 'greater' or 'less'. Can be provided as a
#' single string, or a character vector of the same length as the number of features
#' @param feat.sel.max.pvalue Only features with p-value lower than this will be kept
#' @param feat.sel.min.cliff.delta Cliff delta (effect size) threshold for keeping continous features
#' @param feat.sel.min.cramer.v Cramer's V (effect size) threshold for keeping categorical features
#' @param feat.sel.top.n.features Only keep the top number of features
#' @param balance.classes If 'resample', up and downsampling will be performed. Resampling ratios
#' are determined by an inner cross-validation. If 'class_weights', class weights equal to
#' 1/n_samples_in_class will be used. If 'none', no class balancing will be performed.
#' @param k.inner Number of inner cross-validation folds to use to determine resampling ratios
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
#' @param feat.tmp.path Tmp path for feature selection output
#' @param model.tmp.path Tmp path for trained random forest
#' @param seed Random seed as an integer
#' @param msg.prefix Add prefix to progress messages. If not specified, no prefix will be displayed.
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#' @param feat.sel.whitelist A character vector of feature names in x to keep (i.e. ignore feature
#' selection for these features)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForest <- function(
   train, test=NULL, colname.response='response',

   ## univariate feature slection
   do.feat.sel=T,
   feat.sel.alternative=NULL,
   feat.sel.max.pvalue=0.01,
   feat.sel.min.cliff.delta=0.1,
   feat.sel.min.cramer.v=0.1,
   feat.sel.top.n.features=NULL,
   feat.sel.whitelist=NULL,

   ## class balancing
   balance.classes=c('none','resample','class_weights'), k.inner=5,
   n.resamples.true=4, n.resamples.false=4, midpoint.type='geometric',
   min.size.diff=NULL, max.upsample.ratio=10,

   ntree=500, get.local.increments=T,
   feat.tmp.path=NULL,
   model.tmp.path=NULL,
   seed=NULL, msg.prefix=NULL, verbose=1
){
   if(F){
      #train=folds[[1]]$train
      #test=folds[[1]]$test
      colname.response='response'

      ## univariate feature slection
      do.feat.sel=T; feat.sel.alternative=NULL;
      feat.sel.max.pvalue=0.01;
      feat.sel.min.cliff.delta=0.1;
      feat.sel.min.cramer.v=0.1;
      feat.sel.top.n.features=300;
      feat.sel.whitelist=NULL;

      ## class balancing
      balance.classes=c('none','resample','class_weights'); k.inner=5;
      n.resamples.true=4; n.resamples.false=4; midpoint.type='geometric';
      min.size.diff=NULL; max.upsample.ratio=10;

      ntree=500; get.local.increments=T;
      seed=NULL; msg.prefix=NULL; verbose=3
   }

   if(!is.null(seed)){ set.seed(seed) }

   msg_prefix <- if(!is.null(msg.prefix)){ msg.prefix } else { '' }

   ##----------------------------------------------------------------------
   if(!is.null(model.tmp.path) && file.exists(model.tmp.path)){
      if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Loading saved random forest...') }
      return(readRDS(model.tmp.path))
   }

   ##----------------------------------------------------------------------
   if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Preparing input features and response variable...') }
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

   balance.classes <- balance.classes[1]
   if(!(balance.classes %in% c('none','resample','class_weights'))){
      stop("`balance.classes` must be 'none', 'resample', or 'class_weights'")
   }

   ##----------------------------------------------------------------------
   y <- train_data$y ## y <- train_data$y=='Lung_NSC'
   x <- train_data$x
   rm(train_data)

   out <- list()

   if(do.feat.sel){
      if(!is.null(feat.tmp.path) && !file.exists(feat.tmp.path)){
         if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Performing feature selection...') }
         out$feat_sel <- univarFeatSel(
            x=x, y=y,
            alternative=feat.sel.alternative,
            max.pvalue=feat.sel.max.pvalue,
            min.cliff.delta=feat.sel.min.cliff.delta,
            min.cramer.v=feat.sel.min.cramer.v,
            sel.top.n.features=feat.sel.top.n.features,
            whitelist=feat.sel.whitelist,
            show.avg=F,
            verbose=(verbose>=3)
         )

         saveRDS(out$feat_sel, feat.tmp.path)
      } else {
         if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Loading saved feature selection data...') }
         out$feat_sel <- readRDS(feat.tmp.path)
      }

      keep_features <- out$feat_sel$feature[ out$feat_sel$is_keep_feature==1 ]
      x <- x[,colnames(x) %in% keep_features,drop=F]
   }

   ##----------------------------------------------------------------------
   if(balance.classes=='resample'){
      if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Balancing classes...') }
      resampling_grid <- resamplingGrid(
         a=sum(y==TRUE), b=sum(y==FALSE),
         breaks.a=n.resamples.true, breaks.b=n.resamples.false,
         min.size.diff=min.size.diff, midpoint.type=midpoint.type, max.upsample.ratio=max.upsample.ratio
      )

      train_new <- structure(cbind(y, x), names=c(colname.response, colnames(x)))

      perfs_inner_cv <- lapply(1:nrow(resampling_grid), function(i){
         #i=1
         size_true <- resampling_grid[i,'size_a']
         size_false <- resampling_grid[i,'size_b']

         if(verbose>=2){
            ratio_true <- round(resampling_grid[i,'ratio_a'], 2)
            ratio_false <- round(resampling_grid[i,'ratio_b'], 2)
            message(msg_prefix,'[',format(Sys.time(), "%X"),'] >> TRUE, FALSE: ',ratio_true,'x, ', ratio_false,'x')
         }

         inner_folds <- mltoolkit::createCvTrainTestSets(train_new, stratify.by.col=colname.response, k=k.inner)

         for(j in 1:length(inner_folds)){
            #j=1
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

      x <- train_new[,colnames(train_new)!=colname.response,drop=F]
      y <- train_new[,colname.response]
      rm(train_new)

      out$resampling_grid <- resampling_grid
   }

   if(balance.classes=='class_weights'){
      #y=c(rep(T,10), rep(F,50))
      class_counts <- if(is.logical(y)){
         table(factor(y, c('TRUE','FALSE')))
      } else {
         table(y)
      }
      class_counts <- unclass(class_counts)
      class_weights <- 1/class_counts
   } else {
      class_weights <- NULL
   }

   ##----------------------------------------------------------------------
   if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Training random forest...') }
   model <- randomForest::randomForest(
      x=x,
      y=if(is.logical(y)){ factor(y, c('TRUE','FALSE')) } else { y }, strata=y,
      proximity=F, ntree=ntree, importance=T,
      classwt=class_weights,
      keep.inbag=T, replace=F, ## required for calculating local increments
      na.action=na.roughfix,
      do.trace=F
   )

   if(get.local.increments){
      if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Getting local increments...') }
      ## Used for calculating per sample feature contributions
      invisible(capture.output(
         model$localIncrements <- rfFC::getLocalIncrements(model, x)
      ))
   }

   ## Store output from before training model ------------------------------------------------------
   model$feat_sel <- out$feat_sel
   model$resampling_grid <- out$resampling_grid
   rm(out)

   # ## Store levels from categorical features
   # categorical_lvls <- lapply(x, levels)
   # categorical_lvls <- categorical_lvls[ !sapply(categorical_lvls, is.null) ]
   # categorical_lvls <- categorical_lvls[names(categorical_lvls) %in% colnames(x)]
   #
   # model$categorical_lvls <- categorical_lvls
   # rm(categorical_lvls)

   # ##----------------------------------------------------------------------
   # if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Calculating feature importance...') }
   # model$imp <- sort(
   #    #randomForest::importance(model, type=1)[,1],
   #    model$importance[,'MeanDecreaseAccuracy'],
   #    decreasing=T
   # )

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose>=1){
         message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Predicting on test set...')
      }
      model$test_set <- list(
         probabilities = randomForest:::predict.randomForest(model, test_data$x, type='prob'),
         predicted = NA,
         actual = test_data$y
      )

      model$test_set$predicted <- with(model$test_set,{
         factor(
            colnames(probabilities)[max.col(probabilities)],
            levels=colnames(probabilities)
         )
      })
   }

   model$seed <- seed

   if(!is.null(model.tmp.path)){
      if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Saving tmp model...') }
      saveRDS(model, model.tmp.path)
   }

   return(model)
}

##----------------------------------------------------------------------
#' Train a random forest ensemble (for multiclass classification)
#'
#' @param train A dataframe containing the training features and response column.
#' @param test A dataframe containing the training features and response column. Used for testing
#' the performance of the trained model
#' @param colname.response The column name of the response variable (i.e. training labels)
#' @param do.rmd.nmf Perform NMF to convert RMD bins to RMD signatures?
#' @param args.trainRandomForest A list of arguments that can be passed to `trainRandomForest()`
#' @param seed Random seed as an integer
#' @param fold.num Cross validation fold number. Only used in the prefix of progress messages. If
#' not specified, no prefix will be displayed.
#' @param tmp.dir A path a temporary directory, which if provided, enables resume capability
#' @param rm.tmp.dir If TRUE, will remove the tmp dir if training completes successfully
#' @param multi.core If TRUE, multiple cores will be used when performing NMF on RMD features, and
#' when training the random forest ensemble. Each class (i.e. cancer type) is forked to a separate
#' core.
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForestEnsemble <- function(
   train, test=NULL, colname.response='response',
   do.rmd.nmf=F,
   args.trainRandomForest=list(),
   #inner.holdout.fraction=c(1,3),
   seed=NULL, msg.prefix=NULL,
   tmp.dir=NULL, rm.tmp.dir=T,
   multi.core=F, verbose=1
){
   if(F){
      #set.seed(1)
      #folds=createCvTrainTestSets(training_data)
      #fold=1
      #train=folds[[fold]]$train
      #test=folds[[fold]]$test
      colname.response='response'

      inner.holdout.fraction=c(1,3)
      do.rmd.nmf=T
      #tmp.dir=paste0(dirname(out_path),'/tmp/')
      args.trainRandomForest=list(
         ntree=500,
         do.feat.sel=T,
         feat.sel.alternative=alternative,
         feat.sel.max.pvalue=0.01,
         feat.sel.min.cliff.delta=0.1,
         feat.sel.min.cramer.v=0.1,
         feat.sel.top.n.features=300,
         balance.classes=T,
         verbose=2
      )
      tmp.dir=paste0(dirname(out_path),'/tmp/')

      multi.core=F

      #seed=NULL
      msg.prefix=''
      verbose=2
   }


   ## Init --------------------------------------------------------
   if(!is.null(seed)){ set.seed(seed) }

   msg_prefix <- if(!is.null(msg.prefix)){ msg.prefix } else { '' }

   if(!is.null(tmp.dir)){
      dir.create(tmp.dir, showWarnings=F, recursive=T)
   }

   ## Output list
   out <- list()
   class(out) <- c('list','randomForestEnsemble')
   out$ensemble <- NA
   #out$prob_weigher <- NA

   ##----------------------------------------------------------------------
   if(verbose>=1){message(msg_prefix,'[',format(Sys.time(), "%X"),'] Selecting data for training main ensemble and prob weigher...')}
   #inner.holdout.fraction=c(1,4)

   ## Convert labels to factors
   train[,colname.response] <- as.factor( train[,colname.response] )
   if(any(sapply(train, is.character))){
      stop('Categorical features must be factors and not characters')
   }

   if(!is.null(test)){
      test[,colname.response] <- factor(
         test[,colname.response],
         levels=sort(unique(train[,colname.response]))
      )
      if(any(sapply(test, is.character))){
         stop('Categorical features must be factors and not characters')
      }
   }

   ##----------------------------------------------------------------------
   out$rmd_sig_profiles <- NA

   if(do.rmd.nmf){
      if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Performing NMF on RMD features...') }

      if(!is.null(tmp.dir)){ dir.create(paste0(tmp.dir,'/nmf/'), showWarnings=F) }

      rmd_sig_profiles <- nmfPerClass(
         A=train[,grep('^rmd',colnames(train))],
         response=train[,colname.response],
         tmp.dir=paste0(tmp.dir,'/nmf/'),
         multi.core=multi.core,
         verbose=verbose
      )

      out$rmd_sig_profiles <- t(rmd_sig_profiles$sig_profiles)
      rm(rmd_sig_profiles)

      train <- (function(){
         m1 <- out$rmd_sig_profiles
         m2 <- t(train[,rownames(m1)])

         fit <- NNLM::nnlm(m1, m2)
         fit <- as.data.frame(t(fit$coefficients))
         colnames(fit) <- paste0('rmd.', colnames(fit))

         cbind(
            response=train[,colname.response],
            fit,
            rmColumns(
               train,
               c(grep('^rmd',colnames(train),value=T), colname.response)
            )
         )
      })()

      args.trainRandomForest <- within(args.trainRandomForest,{
         feat.sel.alternative <- feat.sel.alternative[ colnames(train)[-1] ]
         names(feat.sel.alternative) <- colnames(train)[-1]
         feat.sel.alternative[is.na(feat.sel.alternative)] <- 'greater'
      })
   }

   ##----------------------------------------------------------------------
   #y_ohe <- oneHotEncode(as.factor(train_data$ensemble[,colname.response]))
   y_ohe <- oneHotEncode(as.factor(train[,colname.response]))
   #y_ohe <- y_ohe[,c('Lymphoid','Prostate')]

   if(!is.null(tmp.dir)){ dir.create(paste0(tmp.dir,'/class_models/'), showWarnings=F) }
   doTrain <- function(i){
      #i=8
      #i='Gallbladder'

      tmp.feat <- paste0(tmp.dir,'/class_models/features.',colnames(y_ohe)[i],'.rds')
      tmp.class_model <- paste0(tmp.dir,'/class_models/model.',colnames(y_ohe)[i],'.rds')

      y <- unname(y_ohe[,i])
      #train_new <- train_data$ensemble
      train_new <- train
      train_new$response <- y_ohe[,i]

      msg_prefix_ct <- paste0('[',i,'/',ncol(y_ohe),': ',colnames(y_ohe)[i],'] ')
      rf <- do.call(
         trainRandomForest,
         c(
            list(
               train=train_new,
               msg.prefix=msg_prefix_ct,
               model.tmp.path=tmp.class_model,
               feat.tmp.path=tmp.feat
            ),
            args.trainRandomForest
         )
      )

      #if(!is.null(tmp.dir)){ saveRDS(rf, tmp_rds) }

      return(rf)
   }

   tmp.ensemble <- paste0(tmp.dir,'/ensemble.rds')
   if(!file.exists(tmp.ensemble)){

      if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Training random forest ensemble...') }

      if(!multi.core){
         out$ensemble <- lapply(1:ncol(y_ohe), doTrain)
      } else {
         require(foreach)
         require(parallel)
         require(doParallel)

         cores <- parallel::detectCores()
         if(verbose){ message('Multicore: ',cores,' cores') }

         registerDoParallel(cores=cores)
         out$ensemble <- foreach(i=1:ncol(y_ohe)) %dopar% { doTrain(i) }
      }
      names(out$ensemble) <- colnames(y_ohe)

      saveRDS(out$ensemble, tmp.ensemble)
   } else {
      if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Loading random forest ensemble...') }
      out$ensemble <- readRDS(tmp.ensemble)
   }

   #train_data$ensemble <- NULL

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Calculating feature importance...') }
   exist_features <- unique(unlist(lapply(out$ensemble, function(model){
      names(model$forest$ncat)
   })))
   features <- colnames(train)[!(colnames(train) %in% colname.response)]
   exist_features <- features[ features %in% exist_features ] ## Preserve original feature order
   rm(features)

   out$imp <- do.call(rbind, lapply(out$ensemble, function(model){
      #i=out$ensemble[[1]]
      df <- randomForest::importance(model, type=1)
      v <- structure(df[,1],names=rownames(df))
      v <- v[exist_features]
      names(v) <- exist_features
      v[is.na(v)] <- 0
      return(v)
   }))

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Calculating feature stats per class...') }
   out$feat_stats <- featStatsPerClass(
      x=rmColumns(train, colname.response),
      y=train[,colname.response],
      verbose=(verbose==2)
   )

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Predicting on test set...') }
      out$test_set <- list(
         probabilities = predict.randomForestEnsemble(
            out,
            rmColumns(test, colname.response),
            type='prob'
         ),
         predicted = NA,
         actual = test[,colname.response]
      )

      out$test_set$predicted <- with(out$test_set,{
         factor(
            colnames(probabilities)[max.col(probabilities)],
            levels=colnames(probabilities)
         )
      })
   }

   out$seed <- seed

   if(rm.tmp.dir & !is.null(tmp.dir)){
      unlink(tmp.dir, recursive=T)
   }

   return(out)
}

##----------------------------------------------------------------------
#' @export
print.randomForestEnsemble <- function(object){
   cat('$ensemble\nBinary RF names:\n')
   print(names(object$ensemble))

   cat('\nOther objects in list:\n')
   cat( paste0('$',names(object)[3:length(object)]) )
}



