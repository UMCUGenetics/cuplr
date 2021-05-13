#' Apply a summary function to a list of matrixes of the same dimensions
#'
#' @param l A list of matrices
#' @param func A summary function (default is mean())
#' @param as.matrix If TRUE, will return a single matrix with the aggregated values. If FALSE, will
#' return a melted dataframe.
#' @param value.name A custom column name for the value column of the melted dataframe
#'
#' @export
#'
aggregateMatrixList <- function(l, func=function(x){ mean(x, na.rm=T) }, as.matrix=F, value.name=NULL){
   #l=lapply(cv_out,`[[`,'imp')

   if(
      !is.list(l) |
      !all( sapply(l,is.matrix) | sapply(l,is.data.frame) )
   ){
      stop('`l` must be a list of matrices or dataframes')
   }

   all_cols <- unique(unlist(lapply(l, colnames)))
   l_melt <- lapply(l, function(i){
      #i=l[[3]]
      df <- as.data.frame(i)
      missing_cols <- all_cols[!(all_cols %in% colnames(df))]
      df[,missing_cols] <- NA
      df_melt <- reshape2::melt(as.matrix(df))
      colnames(df_melt) <- c('class','feature','value')

      return(df_melt)
   })

   df_pre <- suppressWarnings({
      Reduce(function(x,y){
         merge(x,y,by=c(1,2), all=T, sort=F)
      }, l_melt)
   })

   df <- df_pre[,c(1,2)]
   m_values <- as.matrix(df_pre[,-c(1,2)])
   df$value <- apply(m_values,1,function(i){ func(i) })

   if(as.matrix){
      out <- reshape2::acast(df, class ~ feature)
      out <- out[,colnames(l[[1]])]

      return(out)
   } else {
      if(!is.null(value.name)){ colnames(df)[3] <- value.name }
      return(df)
   }
}
