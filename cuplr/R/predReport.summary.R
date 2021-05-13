#' Summarize prediction report
#'
#' @param report A list with the objects with the names: probs_raw, prob,
#' class_pred, class_actual, feat_contrib, imp
#' @param show.which.top.classes An integer vector specifying the top n predicted classes to show
#' @param show.which.top.features An integer vector specifying the top features contributing to the
#' top predicted class to show
#'
#' @return A dataframe
#' @export
#'
predReportSummary <- function(
   report,
   show.which.top.classes=NULL, show.top.class.probs=T,
   show.which.top.features=NULL, show.top.features.contribs=T
){

   if(F){
      report=pred_reports.cv$HMF_PCAWG
      show.which.top.classes=1:3
      show.which.top.features=1:3
   }

   ## Checks --------------------------------
   if(!is.null(show.which.top.classes) & !is.integer(show.which.top.classes)){
      stop('`show.which.top.classes` must be an integer vector')
   }

   if(!is.null(show.which.top.features) & !is.integer(show.which.top.features)){
      stop('`show.which.top.features` must be an integer vector')
   }

   ## Preds --------------------------------
   df <- data.frame(
      sample=rownames(report$prob),
      row.names=NULL
   )

   if('class_actual' %in% names(report)){
      df$class_actual <- report$class_actual
   }

   df$class_pred <- report$class_pred
   df$class_prob <- matrixStats::rowMaxs(report$prob)

   ## Top predicted class --------------------------------
   if(!is.null(show.which.top.classes)){

      m_next_top_n_classes <- t(apply(report$prob,1,function(i){
         i <- sort(i, decreasing=T)
         i <- round(i, 3)
         i <- i[show.which.top.classes]

         if(show.top.class.probs){
            return( paste0(names(i),'=',i) )
         }
         return( names(i) )
      }))
      colnames(m_next_top_n_classes) <- paste0('class_pred.',show.which.top.classes)

      df <- cbind(df, m_next_top_n_classes)
   }

   ## Top contributing features --------------------------------
   if(!is.null(show.which.top.features)){

      feat_contrib <- report$feat_contrib
      feat_contrib <- feat_contrib[
         paste0(feat_contrib$sample,'_',feat_contrib$binary_rf) %in%
            paste0(df$sample,'_',df$class_pred)
         ,]

      rle_out <- rle(
         paste0(feat_contrib$sample,'_',feat_contrib$binary_rf)
      )
      feat_contrib$index <- unlist(lapply(rle_out$lengths, function(i){ 1:i }))
      feat_contrib <- feat_contrib[
         feat_contrib$index %in% show.which.top.features,
         c('sample','feature','index','contrib')
      ]

      if(show.top.features.contribs){
         feat_contrib$string <- paste0(
            feat_contrib$feature,'=',
            round(feat_contrib$contrib,3)
         )
      } else {
         feat_contrib$string <- feat_contrib$feature
      }

      m_feat_contrib <- cast(feat_contrib, row.var='sample',col.var='index', value.var='string')
      colnames(m_feat_contrib) <- paste0('feat.',colnames(m_feat_contrib))
      m_feat_contrib <- m_feat_contrib[df$sample,]

      df <- cbind(df, m_feat_contrib)
   }

   rownames(df) <- NULL
   return(df)
}
