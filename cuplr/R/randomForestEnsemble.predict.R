#' Predict method for random forest ensemble
#'
#' @param object An object of class randomForestEnsemble
#' @param newdata A data frame or matrix containing new data. (Note: If not given, the out-of-bag
#' prediction in object is returned)
#' @param type 'class', 'prob' or 'votes', indicating the type of output: predicted values,
#' matrix of class probabilities, or matrix of vote counts.
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
   gender.feature.name='gender.gender',
   classes.female=c('Cervix','Ovary','Uterus'), classes.male='Prostate',
   calc.feat.contrib=T, top.n.pred.classes=NULL, top.n.features=5,
   verbose=F
){
   # if(F){
   #    object=model
   #    newdata=features
   #    type='prob'
   #    gender.feature.name='gender.gender'
   #    classes.female=c('Cervix','Ovary','Uterus')
   #    classes.male='Prostate'
   #    top.n.pred.classes=NULL
   #    top.n.features=5
   #    verbose=T
   # }

   ## Checks --------------------------------
   if(!is.data.frame(newdata)){ stop('`newdata` must be a dataframe') }
   if(!('randomForestEnsemble' %in% class(object))){ stop('`object` must be randomForestEnsemble') }

   if(!(type %in% c('prob','class','report'))){
      stop('`type` must be one of the following: prob, class, report')
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
   #probs_adjusted <- randomForest:::predict.randomForest(object$prob_weigher, probs_raw, type='prob')
   probs_adjusted <- probs_raw

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
   classes_pred <- factor( colnames(probs_adjusted)[ max.col(probs_adjusted) ], levels=colnames(probs_adjusted) )
   if(type=='class'){ return(classes_pred) }

   out <- list(
      prob=probs_adjusted,
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

   # m1 <- do.call(cbind, lapply(l_feat_contrib, rowSums))
   # probs_raw

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

   ## --------------------------------
   class(out) <- c('predReport', class(out))
   return(out)
}



