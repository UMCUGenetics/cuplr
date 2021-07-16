#' Predict method for random forest ensemble
#'
#' @param object An object of class randomForestEnsemble
#' @param newdata A data frame or matrix containing new data. (Note: If not given, the out-of-bag
#' prediction in object is returned)
#' @param type 'class', 'prob' or 'votes', indicating the type of output: predicted values,
#' matrix of class probabilities, or matrix of vote counts
#' @param prob.cal.curves A dataframe containing the calibration curve coordinates for scaling the
#' raw probabilities outputted by the classifier (with the column names: x, y, class). Scaled
#' probabilities directly represent accuracy of prediction (i.e. probability of correct prediction)
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
#' @return 'class': predicted classes (the classes with majority vote).
#' 'prob': matrix of class
#' probabilities (one column for each class and one row for each input)
#' 'report': a list containing the following objects: probs_raw, probs_adjusted, responses_pred,
#' feat_contrib
#'
#' @export
#'
predict.randomForestEnsemble <- function(
   object, newdata, type='report',
   prob.cal.curves=NULL,
   gender.feature.name='gender.gender',
   classes.female=c('Cervix','Ovary','Uterus'), classes.male='Prostate',
   calc.feat.contrib=T, top.n.pred.classes=NULL, top.n.features=5,
   verbose=F
){
   if(F){
      object=model
      newdata=features
      type='report'
      prob.cal.curves=cal_curves
      gender.feature.name='gender.gender'
      classes.female=c('Cervix','Ovary','Uterus')
      classes.male='Prostate'
      top.n.pred.classes=NULL
      top.n.features=5
      verbose=T
   }

   ## Checks --------------------------------
   if(!is.data.frame(newdata)){ stop('`newdata` must be a dataframe') }
   if(!('randomForestEnsemble' %in% class(object))){ stop('`object` must be randomForestEnsemble') }

   if(!(type %in% c('prob','class','report'))){
      stop('`type` must be one of the following: prob, class, report')
   }

   if(!is.null(prob.cal.curves)){
      if(!is.data.frame(prob.cal.curves)){
         stop('`prob.cal.curves` must be a dataframe with the columns: x, y, class')
      }
   }

   ## Prepare data --------------------------------
   if(is.matrix(object$rmd_sig_profiles)){
      if(verbose){ message('Fitting RMD profiles...') }
      newdata <- (function(){
         fit <- mutSigExtractor::fitToSignatures(
            mut.context.counts=newdata[,rownames(object$rmd_sig_profiles)], ## Select the columns corresponding to RMD bins
            signature.profiles=object$rmd_sig_profiles
         )
         colnames(fit) <- paste0('rmd.',colnames(fit))
         cbind(
            fit,
            rmColumns( newdata, grep('^rmd',colnames(newdata), value=T) ) ## Remove RMD bin columns
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
   if(verbose){ message('Adjusting raw probabilities based on sample gender...') }

   ## Set probs for disallowed tissue classes to 0
   samples_female <- newdata[,gender.feature.name]
   samples_male <- !samples_female

   probs_adjusted <- probs_raw
   probs_adjusted[samples_female, classes.male] <- 0
   probs_adjusted[samples_male, classes.female] <- 0

   # ## Adjusted probs to sum to 1
   # ## In theory not necessary but leads to better calibrated prob performance
   # probs_adjusted <- probs_adjusted / rowSums(probs_adjusted)

   if(is.null(prob.cal.curves)){
      probs <- probs_adjusted
      prob_scaled <- NULL
   } else {
      if(verbose){ message('Calibrating probabilities...') }
      class(prob.cal.curves) <- c('isoReg',class(prob.cal.curves))
      prob_scaled <- probCal(probs=probs_adjusted, cal.curves=prob.cal.curves)
      probs <- prob_scaled
   }

   if(type=='prob'){ return(probs) }

   ## --------------------------------
   classes_pred <- factor( colnames(probs)[ max.col(probs) ], levels=colnames(probs) )
   if(type=='class'){ return(classes_pred) }

   out <- list(
      prob=probs_adjusted,
      prob_scaled=prob_scaled,
      class_pred=classes_pred
   )

   if(!calc.feat.contrib){ return(out) }

   ## --------------------------------
   if(verbose){
      counter <- 0
      if(verbose){ message('Calculating feature contributions...') }
   }
   l_feat_contrib <- lapply(object$ensemble, function(model){
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

   ## --------------------------------
   if(verbose){ message('Formatting feature contrib output as long form dataframe...') }
   feat_contrib <- reshape2::melt(l_feat_contrib)
   colnames(feat_contrib) <- c('sample','feature','contrib','binary_rf')
   feat_contrib <- feat_contrib[,c('sample','binary_rf','feature','contrib')]

   ## Get features used by ensemble
   all_binary_rf_feat <- unique(unlist(lapply(object$ensemble, function(i){ names(i$forest$xlevels) }), use.names=F))
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

   ## Reduce dataframe size by removing unused data
   if(!is.null(top.n.features)){
      feat_contrib <- feat_contrib[feat_contrib$feature_rank<=top.n.features,]
   }

   if(!is.null(top.n.pred.classes)){
      top_pred_classes <- t(apply(probs,1, function(i){
         colnames(probs)[ order(i, decreasing=T) ]
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

   ## --------------------------------
   if(verbose){ message('Gathering feature summary stats...') }

   ## Get feature values per sample
   newdata_ss <- newdata[,unique(as.character(feat_contrib$feature)),drop=F] ## Subset for features in `feat_contrib`
   newdata_ss <- as.matrix(newdata_ss+0) ## Convert logical data to integer
   newdata_ss <- reshape2::melt(newdata_ss)
   colnames(newdata_ss) <- c('sample','feature','value')

   index <- match(
      paste0(feat_contrib$sample,':',feat_contrib$feature),
      paste0(newdata_ss$sample,':',newdata_ss$feature)
   )
   feat_contrib$sample_value <- newdata_ss$value[index]
   rm(index, newdata_ss)

   ## Get feature statistics per class
   feat_stats <- object$feat_stats
   index <- match(
      paste0(feat_contrib$binary_rf,':',feat_contrib$feature),
      paste0(feat_stats$class,':',feat_stats$feature)
   )

   feat_contrib <- cbind(
      feat_contrib,
      feat_stats[index,c('min_all','max_all','avg_case','avg_ctrl','avg_metric')]
   )
   rownames(feat_contrib) <- NULL

   # index <-
   #    paste0(feat_stats$class,':',feat_stats$feature) %in%
   #    paste0(feat_contrib$binary_rf,':',feat_contrib$feature)
   #
   # feat_stats <- feat_stats[index,]

   out$feat_contrib <- feat_contrib
   #out$feat_stats <- feat_stats
   #rm(feat_stats)

   ## --------------------------------
   class(out) <- c('predReport', class(out))
   return(out)
}

#' @export
print.predReport <- function(object){
   cat('Objects in list:\n')
   cat( paste0('$',names(object)) )

   cat('\n\n')
   if(!is.null(object$prob_scaled)){
      cat('Calibrated probabilities:\n$prob_scaled\n')
      print(object$prob_scaled)
   } else {
      cat('Raw probabilities:\n$prob\n')
      print(object$prob)
   }
}
