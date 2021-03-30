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

   devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/commonUtils/'))
   devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
   library(mltoolkit)

   features <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10f_rmd1mbSigs/features/features.txt.gz')
   features$response <- as.factor(features$response)
   fold_indexes <- readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10f_rmd1mbSigs/cv_out/01/fold_indexes.rds')
}


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
#' @param msg.prefix Add prefix to progress messages. If not specified, no prefix will be displayed.
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
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
   balance.classes=F, k.inner=5,
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
      feat.sel.min.effect.size=0.1;
      feat.sel.top.n.features=300;

      ## class balancing
      balance.classes=F; k.inner=5;
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

   ##----------------------------------------------------------------------
   y <- train_data$y
   x <- train_data$x
   rm(train_data)

   out <- list()

   if(do.feat.sel){
      if(!is.null(feat.tmp.path) && !file.exists(feat.tmp.path)){
         if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Performing feature selection...') }
         out$feat_sel <- univarFeatSel(
            x, y,
            alternative=feat.sel.alternative,
            max.pvalue=feat.sel.max.pvalue,
            min.cliff.delta=feat.sel.min.cliff.delta,
            min.cramer.v=feat.sel.min.cramer.v,
            sel.top.n.features=feat.sel.top.n.features,
            whitelist=feat.sel.whitelist,
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
   if(balance.classes){
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

   ##----------------------------------------------------------------------
   if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Training random forest...') }
   model <- randomForest::randomForest(
      x=x,
      y=if(is.logical(y)){ factor(y, c('TRUE','FALSE')) } else { y }, strata=y,
      proximity=F, ntree=ntree, importance=T,
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

   ## Store levels from categorical features
   categorical_lvls <- lapply(x, levels)
   categorical_lvls <- categorical_lvls[ !sapply(categorical_lvls, is.null) ]
   categorical_lvls <- categorical_lvls[names(categorical_lvls) %in% colnames(x)]

   model$categorical_lvls <- categorical_lvls
   rm(categorical_lvls)

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

   if(verbose>=1){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] > Saving tmp model...') }
   if(!is.null(model.tmp.path)){
      saveRDS(model, model.tmp.path)
   }

   return(model)
}

####################################################################################################
#' Predict method for random forest ensemble
#'
#' @param object An object of class randomForestEnsemble
#' @param newdata A data frame or matrix containing new data. (Note: If not given, the out-of-bag
#' prediction in object is returned)
#' @param type 'response', 'prob' or 'votes', indicating the type of output: predicted values,
#' matrix of class probabilities, or matrix of vote counts. 'class' is allowed, but automatically
#' converted to "response", for backward compatibility.
#' @param gender.feature.name The name of the feature specifying the gender of each sample
#' @param classes.female A character vector of female cancer type classes. For male samples (as
#' determined by `gender.feature.name`), the probabilities outputted by the binary random forests
#' corresponding to `classes.female` will be set to zero
#' @param classes.male A character vector of male cancer type classes. Similar to `classes.female`
#' @param calc.feat.contrib Calculate feature contributions?
#' @param top.n.pred.classes Number of top predicted classes (i.e. binary random forest names) to
#' keep features from when calculating feature contributions
#' @param top.n.features Number of top features for each sample to keep when calculating feature
#' contributions
#' @param verbose Show progress messages? Can be 0,1,2
#'
#'
#' @return 'response': predicted classes (the classes with majority vote).
#' 'prob': matrix of class
#' probabilities (one column for each class and one row for each input)
#' 'detailed': a list containing the following objects: probs_raw, probs_adjusted, responses_pred,
#' feat_contrib
#'
#' @export
#'
predict.randomForestEnsemble <- function(
   object, newdata, type='response',
   gender.feature.name='purple.is_female',
   classes.female=c('Cervix','Ovary','Uterus'), classes.male='Prostate',
   calc.feat.contrib=T, top.n.pred.classes=3, top.n.features=5,
   verbose=F
){
   if(F){
      object=model
      newdata=features
      #newdata=unlist(features[1,])
      type='prob'
      gender.feature.name='purple.is_female'
      classes.female=c('Cervix','Ovary','Uterus')
      classes.male='Prostate'
      top.n.pred.classes=3
      top.n.features=5
      verbose=T
   }

   ## Checks --------------------------------
   if(!is.data.frame(newdata)){ stop('`newdata` must be a dataframe') }
   if(!('randomForestEnsemble' %in% class(object))){ stop('`object` must be randomForestEnsemble') }

   if(!(type %in% c('prob','report','detailed'))){
      stop('`type` must be one of the following: prob, report, detailed')
   }

   ## Prepare data --------------------------------
   categorical_lvls <- unname(lapply(object$ensemble, function(i){ i$categorical_lvls }))
   categorical_lvls <- unlist(categorical_lvls, recursive=F)
   categorical_lvls <- categorical_lvls[!duplicated(categorical_lvls)]

   if(length(categorical_lvls)>1){
      if(verbose){ message('Assigning categorical variable levels...') }
      for(i in names(categorical_lvls)){
         #i='purple.gender'
         newdata[,i] <- factor(newdata[,i], levels=categorical_lvls[[i]])
      }
   }

   if(is.matrix(object$rmd_sig_profiles)){
      if(verbose){ message('Fitting RMD profiles...') }
      newdata <- (function(){
         m1 <- object$rmd_sig_profiles   ## rows: bins, cols: signatures
         m2 <- t(newdata[,rownames(m1)]) ## rows: bins, cols: samples

         fit <- NNLM::nnlm(m1, m2)
         fit <- as.data.frame(t(fit$coefficients))
         colnames(fit) <- paste0('rmd.',colnames(fit))
         cbind(
            fit,
            rmColumns( newdata, grep('^rmd',colnames(newdata), value=T) )
         )
      })()
   }

   ## --------------------------------
   if(verbose){
      counter <- 0
      message('Getting predictions from each binary RF...')
   }
   #probs_raw <- reports_merged$probs_raw
   probs_raw <- do.call(cbind, lapply(object$ensemble, function(model){
      if(verbose>=2){
         counter <<- counter + 1
         message('[',counter,'/',length(object$ensemble),']: ', names(object$ensemble)[counter])
      }
      randomForest:::predict.randomForest(model, newdata, type='prob')[,1]
   }))

   ## --------------------------------
   if(verbose){ message('Adjusting raw probabilities...') }

   ## Re-weigh probs
   probs_adjusted <- randomForest:::predict.randomForest(object$prob_weigher, probs_raw, type='prob')

   ## Set probs for disallowed tissue classes to 0
   samples_female <- newdata[,gender.feature.name]
   samples_male <- !samples_female

   probs_adjusted[samples_female, classes.male] <- 0
   probs_adjusted[samples_male, classes.female] <- 0

   ## Rescale probs to sum to 1
   probs_adjusted <- probs_adjusted / rowSums(probs_adjusted)

   decimal_places <- max(nchar(sub('^\\d*[.]','',probs_raw)))
   probs_adjusted <- round(probs_adjusted, decimal_places)

   #summary(probs_adjusted[,'Uterus'])
   if(type=='prob'){ return(probs_adjusted) }

   ## --------------------------------
   responses_pred <- factor( colnames(probs_adjusted)[ max.col(probs_adjusted) ], levels=colnames(probs_adjusted) )
   if(type=='response'){ return(responses_pred) }

   out <- list(
      probs_raw=probs_raw,
      probs_adjusted=probs_adjusted,
      responses_pred=structure(responses_pred, names=rownames(newdata))
   )

   if(!calc.feat.contrib){ return(out) }

   ## --------------------------------
   if(verbose){
      counter <- 0
      if(verbose){ message('Calculating feature contributions...') }
   }
   l_feat_contrib <- lapply(model$ensemble, function(model){
      if(verbose>=2){
         counter <<- counter + 1
         message('[',counter,'/',length(object$ensemble),']: ', names(object$ensemble)[counter])
      }

      rfFC::featureContributions(
         object=model,
         lInc=model$localIncrements,
         dataT=newdata
      )
   })

   # m1 <- do.call(cbind, lapply(l_feat_contrib, rowSums))
   # probs_raw

   ## --------------------------------
   if(verbose){ message('Formatting feature contrib output as long form dataframe...') }
   feat_contrib <- reshape2::melt(l_feat_contrib)
   colnames(feat_contrib) <- c('sample','feature','contrib','binary_rf')
   feat_contrib <- feat_contrib[,c('sample','binary_rf','feature','contrib')]

   ## Get features used by ensemble
   all_binary_rf_feat <- unique(unlist(lapply(model$ensemble, function(i){ names(i$forest$xlevels) }), use.names=F))
   all_binary_rf_feat <- all_binary_rf_feat[ na.exclude(match(colnames(newdata), all_binary_rf_feat)) ]

   ## Convert metadata to factors
   feat_contrib <- within(feat_contrib,{
      sample <- factor(sample, rownames(newdata))
      feature <- factor(feature, all_binary_rf_feat)
      binary_rf <- factor(binary_rf, names(l_feat_contrib))
   })

   ## Calculate feature rank
   feat_contrib <- feat_contrib[
      order(feat_contrib$sample, feat_contrib$binary_rf, -feat_contrib$contrib)
      ,]

   feat_contrib$group <- as.numeric(
      paste0(
         as.integer(feat_contrib$sample),'.',
         as.integer(feat_contrib$binary_rf)
      )
   )

   rle_out <- rle(feat_contrib$group)
   feat_contrib$feature_rank <- unlist(lapply(rle_out$lengths, function(i){ 1:i }))
   feat_contrib$group <- NULL; rm(rle_out)

   ## Reduce object size by removing irrelevant data
   if(!is.null(top.n.features)){
      feat_contrib <- feat_contrib[feat_contrib$feature_rank<=top.n.features,]
   }

   if(!is.null(top.n.pred.classes)){
      top_pred_classes <- t(apply(probs_adjusted,1, function(i){
         colnames(probs_adjusted)[ order(i, decreasing=T) ]
      }))
      top_pred_classes <- top_pred_classes[,1:top.n.pred.classes,drop=F]
      top_pred_classes <- reshape2::melt(top_pred_classes)
      colnames(top_pred_classes) <- c('sample','pred_class_rank','pred_class')

      feat_contrib <- feat_contrib[
         paste0(feat_contrib$binary_rf,':',feat_contrib$sample) %in%
            paste0(top_pred_classes$pred_class,':',top_pred_classes$sample)
         ,]

      rm(top_pred_classes)
   }

   rownames(feat_contrib) <- NULL

   ## --------------------------------
   if(verbose){ message('Gathering feature summary stats...') }
   feat_sel <- do.call(rbind, lapply(names(object$ensemble), function(i){
      #i='Prostate'
      df <- object$ensemble[[i]]$feat_sel
      cbind(binary_rf=i, df)
   }))

   ## Get cohort averages
   index <- match(
      paste0(feat_contrib$binary_rf,':',feat_contrib$feature),
      paste0(feat_sel$binary_rf,':',feat_sel$feature)
   )
   feat_contrib$avg_case <- feat_sel$avg_case[index]
   feat_contrib$avg_ctrl <- feat_sel$avg_ctrl[index]
   rm(index)

   ## Get feature values per sample
   newdata_ss <- newdata[,unique(as.character(feat_contrib$feature)),drop=F] ## Subset for features in `feat_contrib`
   newdata_ss <- as.matrix(newdata_ss+0) ## Convert logical data to integer
   newdata_ss <- reshape2::melt(newdata_ss)
   colnames(newdata_ss) <- c('sample','feature','value')

   index <- match(
      paste0(feat_contrib$sample,':',feat_contrib$feature),
      paste0(newdata_ss$sample,':',newdata_ss$feature)
   )
   feat_contrib$value <- newdata_ss$value[index]
   rm(index, newdata_ss)

   ##
   out$feat_contrib <- feat_contrib
   return(out)
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
#' @param rm.tmp.dir If TRUE will remove the tmp dir if training completes successfully
#' @param multi.core Use multiple cores?
#' @param verbose Show progress? Can be 0, 1, 2 (increasing verbosity)
#'
#' @return A list containing the training output
#' @export
#'
trainRandomForestEnsemble <- function(
   train, test=NULL, colname.response='response',
   do.rmd.nmf=F,
   args.trainRandomForest=list(),
   inner.holdout.fraction=c(1,3),
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
   out$prob_weigher <- NA

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
   if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Splitting data for training main ensemble and probability weigher...') }
   inner_folds_indexes <- mltoolkit::createCvTrainTestSets(
      train, k=inner.holdout.fraction[2],
      stratify.by.col=colname.response, return.data=F
   )
   inner_folds_indexes <- lapply(inner_folds_indexes,`[[`,'test')

   train_data <- list()

   train_data$ensemble <- train[
      sort(unlist(inner_folds_indexes[
         (inner.holdout.fraction[1]+1) : inner.holdout.fraction[2]
      ]))
   ,]

   train_data$prob_weigher <- train[
      sort(unlist(inner_folds_indexes[
         1:inner.holdout.fraction[1]
      ]))
   ,]

   ##----------------------------------------------------------------------
   y_ohe <- oneHotEncode(as.factor(train_data$ensemble[,colname.response]))
   #y_ohe <- y_ohe[,c('Lymphoid','Prostate')]

   if(!is.null(tmp.dir)){ dir.create(paste0(tmp.dir,'/class_models/'), showWarnings=F) }
   doTrain <- function(i){
      #i=8
      #i='Gallbladder'

      tmp.feat <- paste0(tmp.dir,'/class_models/features.',colnames(y_ohe)[i],'.rds')
      tmp.class_model <- paste0(tmp.dir,'/class_models/model.',colnames(y_ohe)[i],'.rds')

      y <- unname(y_ohe[,i])
      train_new <- train_data$ensemble
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

   train_data$ensemble <- NULL

   ##----------------------------------------------------------------------
   if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Training random forest probability weigher...') }

   tmp.prob_weigher <- paste0(tmp.dir,'/prob_weigher.rds')
   if(!file.exists(tmp.prob_weigher)){
      holdout_probs <- do.call(cbind, lapply(out$ensemble, function(model){
         #i=ensemble$Breast
         randomForest:::predict.randomForest(
            model,
            rmColumns(train_data$prob_weigher, colname.response),
            type='prob'
         )[,1]
      }))

      holdout_probs <- cbind(
         response=train_data$prob_weigher[,colname.response],
         as.data.frame(holdout_probs)
      )

      out$prob_weigher <- randomForest::randomForest(
         x=rmColumns(holdout_probs, colname.response),
         y=holdout_probs[,colname.response],
         strata=y,
         proximity=F, ntree=500, importance=T,
         keep.inbag=T, replace=F, ## required for calculating local increments
         na.action=na.roughfix,
         do.trace=F
      )

      invisible(capture.output(
         out$prob_weigher$localIncrements <- rfFC::getLocalIncrements(
            out$prob_weigher,
            rmColumns(holdout_probs, colname.response)
         )
      ))
      saveRDS(out$prob_weigher, tmp.prob_weigher)
      rm(holdout_probs)

   } else {

      if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Loading random forest probability weigher...') }
      out$prob_weigher <- readRDS(tmp.prob_weigher)
   }

   train_data$prob_weigher <- NULL

   # ##----------------------------------------------------------------------
   # if(verbose){ message(msg_prefix,'[',format(Sys.time(), "%X"),'] Storing levels from categorical features...') }
   # categorical_lvls <- lapply(rmColumns(train_data$prob_weigher, colname.response), levels)
   # categorical_lvls <- categorical_lvls[ !sapply(categorical_lvls, is.null) ]
   # categorical_lvls <- categorical_lvls[
   #    names(categorical_lvls) %in% colnames(rmColumns(train_data$prob_weigher, colname.response))
   # ]
   #
   # out$categorical_lvls <- categorical_lvls

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

   ## --------
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
   cat('$ensemble\n\nBinary RF names:\n')
   print(names(object$ensemble))

   cat('\n$prob_weigher')
   print(object$prob_weigher)

   cat('\nOther list levels:\n')
   cat( paste0('$',names(object)[3:length(object)]) )
}



