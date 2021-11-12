#' Summarize prediction report
#'
#' @param report A list with the objects with the names: probs_raw, prob,
#' class_pred, class_actual, feat_contrib, imp
#' @param top.n.classes An integer specifying the top n predicted classes to show
#' @param show.class.probs When `top.n.classes`==TRUE, also show the predicted
#' probabilities?
#' @param simplify.pred.values If TRUE, the top n predicted classes and probs will be combined as a
#' string in the form {class}={prob}
#' @param prob.type Can be 'prob' (raw probabilities) or 'prob_scaled' (scaled probabilities)
#' @param top.n.feat An integer specifying the top features contributing to the top predicted
#' class to show
#' @param top.n.feat.classes An integer specifying the number of top classes for which to show the
#' top contributing features
#' @param show.feat.contribs When `top.n.feat`==TRUE, also show the feature
#' contributions?
#'
#' @return A dataframe
#'
#' @method summary predReport
#' @export
#'
summary.predReport <- function(
   report,
   top.n.classes=NULL, show.class.probs=T, simplify.pred.values=F, prob.type=NULL,
   top.n.feat=NULL, top.n.feat.classes=1, show.feat.contribs=T
){

   if(F){
      report=pred_reports$CUP
      report=pred_report
      top.n.classes=3
      show.class.probs=T
      simplify.pred.values=T
      prob.type=NULL
      top.n.feat.classes=1
      top.n.feat=3
      show.feat.contribs=T
   }

   ## Init --------------------------------
   ##
   if(!is.null(top.n.classes) & !is.numeric(top.n.classes)){
      stop('`top.n.classes` must be an integer vector')
   }

   if(!is.null(top.n.feat) & !is.numeric(top.n.feat)){
      stop('`top.n.feat` must be an integer vector')
   }

   if(is.null(prob.type)){
      if(!is.null(report$prob_scaled)){
         prob.type <- 'prob_scaled'
      } else {
         prob.type <- 'prob'
      }
   }
   report$prob <- report[[prob.type]]

   if(length(rownames(report$prob))==0){
      rownames(report$prob) <- 1:nrow(report$prob)
   }
   df <- data.frame(sample=rownames(report$prob), row.names=NULL)

   ## Top predicted class --------------------------------
   top_classes <- t(apply(report$prob,1,function(i){ names(sort(i, decreasing=T)) }))
   top_probs <- t(apply(report$prob,1,function(i){ sort(i, decreasing=T) }))

   if('class_actual' %in% names(report)){
      df$actual_class <- report$class_actual

      top_correct <- sweep(top_classes, 1, report$class_actual, '==')

      ## Get the top pred class that is the correct class
      df$pred_correct <- max.col(top_correct, ties.method='last')

      ## Assign 0 when actual class doesn't match any of the predicted classes
      df$pred_correct[ rowSums(top_correct)==0 ] <- 0
      #rm(top_correct)
   }

   if(is.null(top.n.classes)){
      df$pred_class.1 <- top_classes[,1]
      df$pred_prob.1 <- top_probs[,1]

   } else {

      if(length(top.n.classes)==1){ top.n.classes <- 1:top.n.classes }

      top_probs <- round(top_probs, 3)
      top_probs <- top_probs[,top.n.classes,drop=F]
      colnames(top_probs) <- paste0('pred_prob.',top.n.classes)

      top_classes <- top_classes[,top.n.classes,drop=F]
      colnames(top_classes) <- paste0('pred_class.',top.n.classes)

      ##
      if(simplify.pred.values){
         pred <- matrix(
            paste0(top_classes,'=',top_probs),
            nrow=nrow(top_probs), ncol=ncol(top_probs),
            dimnames=list( NULL, paste0('pred_class.',top.n.classes) )
         )
      } else {
         pred <- data.frame(top_classes, top_probs, row.names=NULL)
      }

      df <- data.frame(df, pred)
   }
   #rm(top_probs, top_classes)

   ## Top contributing features --------------------------------
   if(!is.null(top.n.feat)){

      if(length(top.n.feat)==1){
         top.n.feat <- 1:top.n.feat
      }

      feat_contrib <- subset(report$feat_contrib, feature_rank %in% top.n.feat)

      m_feat_contrib <- lapply(1:top.n.feat.classes, function(i){
         #i=1

         feat_contrib_ss <- feat_contrib[
            paste0(feat_contrib$sample,'_',feat_contrib$binary_rf) %in%
               paste0(df$sample,'_',top_classes[,i])
            ,]

         ## Get feature rank within each sample
         rle_out <- rle( paste0(feat_contrib_ss$sample,'_',feat_contrib_ss$binary_rf) )
         feat_contrib_ss$index <- unlist(lapply(rle_out$lengths, function(i){ 1:i }))
         feat_contrib_ss <- feat_contrib_ss[
            feat_contrib_ss$index %in% top.n.feat,
            c('sample','feature','index','contrib')
         ]

         ## Make final string output
         if(show.feat.contribs){
            feat_contrib_ss$string <- paste0( feat_contrib_ss$feature, '=', round(feat_contrib_ss$contrib,3) )
         } else {
            feat_contrib_ss$string <- feat_contrib_ss$feature
         }

         ## Convert long to wide form dataframe
         m <- reshape2::dcast(data=feat_contrib_ss, formula=sample~index, value.var='string', fill=NA_real_)
         rownames(m) <- m[,1]; m <- m[,-1]
         m <- m[as.character(df$sample),]

         ## Add column names
         colnames(m) <- paste0('feat_contrib.',i,'.',colnames(m))

         return(m)
      })
      m_feat_contrib <- do.call(cbind, m_feat_contrib)

      df <- data.frame(df, m_feat_contrib, row.names=NULL)
   }

   return(df)
}
