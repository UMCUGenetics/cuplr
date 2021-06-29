#' Summary stats on a dataframe of features for each class
#'
#' @description This function calculates the average feature values between the target cancer type
#' samples vs other samples (interquartile mean for numeric data, or proportion for logical data),
#' as well as min and max feature values over all samples. This data is used for plotting feature
#' contributions in the patient report.
#'
#' @param x A dataframe of features where rows are samples
#' @param y A character or factor vector of class labels for each sample in `x`
#' @param verbose Show progress messages?
#'
#' @return A long form dataframe with the min, max, and average feature values for each feature for
#' each class
#' @export
#'
featStatsPerClass <- function(x, y, verbose=F){

   # if(F){
   #    x=features[,-1]
   #    y=features[,1]
   # }

   require(matrixStats)

   ## Checks --------------------------------
   if(!is.data.frame(x)){ stop('x must be a dataframe') }
   #if(!is.logical(y) | !is.character(y)){ stop('y must be a logical vector') }

   is_numeric <- sapply(x, is.numeric, USE.NAMES=F)
   is_logical <- sapply(x, is.logical, USE.NAMES=F)
   if(!all(is_numeric | is_logical)){ stop('x must only contain numeric or logical data') }

   if(is.null(colnames(x))){ stop('x must have colnames') }

   ## Main function --------------------------------
   if(verbose){ message('Splitting numeric and logical features') }
   x_numeric <- x[,is_numeric,drop=F]
   x_numeric <- as.matrix(x_numeric)

   x_logical <- x[,is_logical,drop=F]
   x_logical <- as.matrix(x_logical)+0 ## Use +0 to convert to integer matrix

   main <- function(x.numeric, x.logical, calc.avg=T, calc.min.max=T){
      #x.numeric=x_numeric
      #x.logical=x_logical

      ## Numeric data
      min_numeric <- numeric()
      max_numeric <- numeric()
      avg_numeric <- numeric()
      avg_metric_numeric <- character()

      if(ncol(x.numeric)!=0){
         if(calc.min.max){
            min_numeric <- colMins(x.numeric)
            max_numeric <- colMaxs(x.numeric)
         }

         if(calc.avg){
            avg_numeric <- colMeansTrimmed(x.numeric, trim=0.25, na.rm=T)
            avg_metric_numeric <- rep('iqm',ncol(x.numeric))
         }
      }

      ## Logical data
      min_logical <- numeric()
      max_logical <- numeric()
      avg_logical <- numeric()
      avg_metric_logical <- character()

      if(ncol(x.logical)!=0){
         if(calc.min.max){
            min_logical <- rep(0,ncol(x.logical))
            max_logical <- rep(1,ncol(x.logical))
         }

         if(calc.avg){
            avg_logical <- colMeans(x.logical, na.rm=T) ## Using colMeans() is equal taking the proportion of TRUE
            avg_metric_logical <- rep('prop',ncol(x.logical))
         }
      }

      ## Gather stats
      out <- data.frame(
         feature=c(colnames(x.numeric), colnames(x.logical))
      )

      if(calc.min.max){
         out$min <- c(min_numeric, min_logical)
         out$max <- c(max_numeric, max_logical)
      }

      if(calc.avg){
         out$avg <- c(avg_numeric, avg_logical)
         out$avg_metric <- c(avg_metric_numeric, avg_metric_logical)
      }

      ## Restore original feature order
      rownames(out) <- out$feature
      out <- out[colnames(x),]
      rownames(out) <- NULL

      return(out)
   }

   ## Apply main function to all samples and each class --------------------------------
   if(verbose){ message('Calculating stats for:') }
   classes <- sort(unique(y))
   df_all <- main(x_numeric, x_logical, calc.avg=F)
   feat_stats <- do.call(rbind, lapply(classes, function(i){
      message('> ',i)
      #i='Prostate'

      df_case <- main(
         x_numeric[y==i,,drop=F],
         x_logical[y==i,,drop=F],
         calc.min.max=F
      )

      df_ctrl <- main(
         x_numeric[y!=i,,drop=F],
         x_logical[y!=i,,drop=F],
         calc.min.max=F
      )

      data.frame(
         class=i,
         feature=df_case$feature,
         min_all=df_all$min,
         max_all=df_all$max,
         avg_case=df_case$avg,
         avg_ctrl=df_ctrl$avg,
         avg_metric=df_case$avg_metric
      )
   }))
}
