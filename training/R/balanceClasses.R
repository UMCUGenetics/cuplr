################################################################################
#' Create pairs of resample sizes/ratios for performing a grid search
#'
#' @param a The sample size of cohort 1
#' @param b The sample size of cohort 1
#' @param breaks.a Number of resampling values to generate for a
#' @param breaks.b Number of resampling values to generate for b
#' @param min.size.diff Default=30. Minimum difference between resampling values (integer).
#' @param midpoint.type Can be 'geometric', 'arithmetic', or 'none'. Calculate the resampling values
#' from a->midpoint and b->midpoint? If 'none', resampling values will be calculated from a->b and
#' b->a.
#' @param max.upsample.ratio Default=10. Remove pairs where a or b are upsampled more than this
#' value
#'
#' @return A dataframe of target sample sizes and resampling ratios for a and b
#' @export
#'
resamplingGrid <- function(
   a, b,
   breaks.a=4, breaks.b=4,
   midpoint.type='geometric',
   min.size.diff=30, max.upsample.ratio=10
){
   # if(F){
   #    a=25
   #    b=3500
   #    breaks.a=5
   #    breaks.b=5
   #    midpoint.type='geometric'
   #    min.size.diff=30
   #
   # }

   ## Calculate midpoint between a and b
   if(midpoint.type=='geometric'){
      midpoint <- round(sqrt(a*b)) ## logarithmic (geometric) mean Same as exp((log(a)+log(b))/2)
   } else if(midpoint.type=='arithmetic') {
      midpoint <- (a + b)/2
   } else {
      midpoint <- NA ## Use a and b as endpoints
   }

   if(is.null(min.size.diff)){
      if(!is.na(midpoint)){
         ## Automatically calculate minimum difference between intervals
         min.size.diff <- min(
            round(sqrt(abs(midpoint-a))),
            round(sqrt(abs(midpoint-b)))
         )
      } else {
         min.size.diff <- 0
      }
   }

   ## Main
   calcTargetSampleSize <- function(start, end, breaks){
      #start=500
      #end=midpoint
      #breaks=4

      if(is.null(min.size.diff) || min.size.diff<=0){
         return(
            round( 2 ^ seq(log2(start), log2(end), length.out=breaks) )
         )
      }

      difference <- 0
      while(difference < min.size.diff){
         size <- 2 ^ seq(log2(start), log2(end), length.out=breaks)
         size <- round(size)
         difference <- min(abs(diff(size)))
         breaks <- breaks - 1
         if(breaks<=1){ break }
      }
      return(size)
   }

   size_a <- calcTargetSampleSize(
      start=a,
      end=if(!is.na(midpoint)){ midpoint } else { b },
      breaks=breaks.a
   )

   size_b <- calcTargetSampleSize(
      start=b,
      end=if(!is.na(midpoint)){ midpoint } else { a },
      breaks=breaks.b
   )

   ## Return output
   out <- structure(
      expand.grid(size_a, size_b),
      names=c('size_a','size_b')
   )
   out$ratio_a <- out$size_a / a
   out$ratio_b <- out$size_b / b

   if(!is.null(max.upsample.ratio)){
      out <- out[
         out$ratio_a < max.upsample.ratio &
         out$ratio_b < max.upsample.ratio
      ,]

   }

   return(out)
}

#resamplingGrid(460, 3500, breaks.a=4, breaks.b=4, min.size.diff=NULL)

################################################################################
#' Resolve class imbalances by simple resampling
#'
#' @param df A dataframe containing the features and response variable
#' @param colname.response Column name of the response variable
#' @param resample.ratios A named numeric vector indicating the ratio of upsampling for each class.
#' Values < 1 indicate downsampling. Values == 1 indicates no resampling
#' @param target.sample.sizes A named integer vector indicating the target number of samples for
#' each class
#' @param return.data If TRUE, returns the resampled dataframe. If FALSE, returns the row indexes.
#'
#' @return A resampled dataframe or a vector of row indexes
#' @export
#'
resampleClasses <- function(
   df, colname.response='response',
   resample.ratios=NULL, target.sample.sizes=NULL,
   return.data=T
){
   #df=train[,1:10]
   #resample.ratios=c('TRUE'=1.5,'FALSE'=1)

   y <- df[,colname.response]
   indexes <- split(1:nrow(df), y)

   if(!is.null(resample.ratios)){
      if(is.null(names(resample.ratios))){
         stop('`resample.ratios` must be a numeric vector with names (of the sample classes)')
      }
      if(!all(names(resample.ratios) %in% names(indexes))){
         stop('Not all class names are present in the names of `resample.ratios`')
      }

      target.sample.sizes <- sapply(names(indexes), function(i){
         round( length(indexes[[i]]) * resample.ratios[[i]] )
      })
   }

   if(is.null(names(target.sample.sizes))){
      stop('`target.sample.sizes` must be an integer vector with names (of the sample classes)')
   }
   if(!all(names(target.sample.sizes) %in% names(indexes))){
      stop('Not all class names are present in the names of `target.sample.sizes`')
   }

   indexes_new <- lapply(names(target.sample.sizes), function(i){
      v <- indexes[[i]]
      target_sample_size <- target.sample.sizes[[i]]

      if(length(v)==target_sample_size){ return(v) }

      sample(
         v,
         size=target_sample_size,
         replace=if(target_sample_size > length(v)){ TRUE } else { FALSE }
      )
   })
   #names(indexes_new) <- names(target.sample.sizes)

   indexes_new <- sort(unlist(indexes_new))

   if(return.data){ df[indexes_new,] } else { indexes_new }
}
