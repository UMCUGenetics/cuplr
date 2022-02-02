################################################################################
#' Create pairs of resample sizes/ratios for performing a grid search
#'
#' @description This function aims to automatically select the best pairs of resampling parameters
#' given the sample size of to groups `a` and `b`. The priority for selection of pairs is as follows:
#' * downsampling and no upsampling
#' * pass `min.downsample.size`
#' * pass `min.downsample.ratio`
#' * pass `max.upsample.size`
#' * pass `max.upsample.ratio`
#' * pass `max.imbalance.ratio`: prioritize pairs resulting in acceptably low imbalances
#' * resampling intensity: prioritizes lower intensities; intensity = (1/ratio_b + 1) ^ ratio_a
#' * size_a: prioritizes smaller a sample sizes
#' * size_b: prioritizes larger b sample sizes
#'
#' @param a The sample size of cohort a
#' @param b The sample size of cohort b. b must be > a
#' @param breaks.a Number of resampling values to generate for a
#' @param breaks.b Number of resampling values to generate for b
#' @param midpoint.type Can be 'geometric', 'arithmetic', or 'none'. Calculate the resampling values
#' from a->midpoint and b->midpoint? If 'none', resampling values will be calculated from a->b and
#' b->a.
#' @param min.size.diff Default=15. Minimum difference between resampling values (integer).
#' @param min.downsample.size Default=1000. Pairs where downsampling of `b` leads to a sample size less than this value are deprioritized
#' @param min.downsample.ratio Default=NULL. Pairs where downsampling of `b` leads to a resampling ratio less than this value are deprioritized
#' @param max.upsample.size Default=300. Pairs where upsampling of `a` leads to a sample size greater than this value are deprioritized
#' @param max.upsample.ratio Default=5. Pairs where upsampling of `a` leads to a resampling ratio greater than this value are deprioritized
#' @param max.imbalance.ratio Default=15. Pairs where resampling leads high `max(a/b, b/a)` are deprioritized
#' @param max.pairs Default=10. The max number of target sample size pairs to return.
#'
#' @param bypass.filter If TRUE, all pairs of target sample sizes are returned. If FALSE, only the
#' pairs strictly passing the follow criteria are returned:
#' * pass `max.upsample.size`
#' * pass `min.downsample.size`
#' * max.pairs
#'
#' @return A dataframe of target sample sizes and resampling ratios for a and b
#' @export
#'
resamplingGrid <- function(
   a, b,
   breaks.a=5, breaks.b=5,
   midpoint.type='geometric',

   min.size.diff=15,
   min.downsample.size=1000, min.downsample.ratio=NULL,
   max.upsample.size=300, max.upsample.ratio=5,
   max.imbalance.ratio=15, max.pairs=10, bypass.filter=F
){
   # if(F){
   #    a=50
   #    b=6000
   #    breaks.a=5
   #    breaks.b=5
   #    midpoint.type='geometric'
   #    max.upsample.ratio=5
   #    max.upsample.size=300
   #    min.downsample.ratio=NULL
   #    min.downsample.size=1000
   #    max.imbalance.ratio=15
   #    max.pairs=10
   # }

   if(a>b){ stop('b must be greater than a') }
   #if(a==b){ }

   ## Calculate midpoint between a and b --------------------------------
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

   ## Main --------------------------------
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

   if(breaks.a==1){
      size_a <- a
   } else {
      size_a <- calcTargetSampleSize(
         start=a,
         end=if(!is.na(midpoint)){ midpoint } else { b },
         breaks=breaks.a
      )
   }

   if(breaks.b==1){
      size_b <- b
   } else {
      size_b <- calcTargetSampleSize(
         start=b,
         end=if(!is.na(midpoint)){ midpoint } else { a },
         breaks=breaks.b
      )
   }

   ## Filtering --------------------------------
   out <- structure(
      expand.grid(size_a, size_b),
      names=c('size_a','size_b')
   )

   ##
   out$size_b[ out$size_b < min.downsample.size ] <- min.downsample.size
   out <- unique(out)

   ##
   out$ratio_a <- out$size_a / a
   out$ratio_b <- out$size_b / b

   out$imbalance_ratio <- pmax(
      out$size_a/out$size_b,
      out$size_b/out$size_a
   )

   ## Calculate resampling intensity score
   ## Use ratio_a as exponent to penalize more upsampling
   ## Add 1 to 1/ratio_b to prevent 1^n (i.e. penalty still applies when no downsampling is performed)
   out$resampling_intensity <- (1/out$ratio_b + 1) ^ out$ratio_a

   ##
   out$pass.min.downsample.size <- TRUE
   if(!is.null(min.downsample.size)){
      ## Never allow too much down sampling
      out$pass.min.downsample.size[out$ratio_b<1 & out$size_b<min.downsample.size] <- FALSE
   }

   out$pass.min.downsample.ratio <- TRUE
   if(!is.null(min.downsample.ratio)){
      out$pass.min.downsample.ratio[out$ratio_b < min.downsample.ratio] <- FALSE
   }

   ##
   out$pass.max.upsample.size <- TRUE
   if(!is.null(max.upsample.size)){
      ## When original sample size is already high, no need to upsample
      out$pass.max.upsample.size[a>max.upsample.size & out$ratio_a>1 & out$size_a>max.upsample.size] <- FALSE
   }

   out$pass.max.upsample.ratio <- TRUE
   if(!is.null(max.upsample.ratio)){
      out$pass.max.upsample.ratio[out$ratio_a > max.upsample.ratio] <- FALSE
   }

   ##
   out$pass.max.imbalance.ratio <- TRUE
   if(!is.null(max.imbalance.ratio)){
      out$pass.max.imbalance.ratio[ out$imbalance_ratio>max.imbalance.ratio ] <- FALSE
   }

   ## Prioritize pairs
   out <- out[
      with(out, order(
         -(out$ratio_a==1 & out$ratio_b<1), ## prioritize downsampling and no upsampling
         -pass.min.downsample.size,
         -pass.min.downsample.ratio,
         -pass.max.upsample.size,
         -pass.max.upsample.ratio,
         -pass.max.imbalance.ratio, ## prioritize pairs resulting in acceptably low imbalances
         resampling_intensity,
         size_a,
         -size_b
      ))
   ,]

   if(!bypass.filter){
      out <- out[out$pass.max.upsample.size,] ## When original sample size is already high, no need to upsample
      out <- out[out$pass.min.downsample.size,] ## Never allow too much down sampling
      out <- out[1:min(max.pairs,nrow(out)),] ## Select maximum number of pairs
      out <- out[,!grepl('^pass',colnames(out))] ## Remove temporary pass columns from output
   }

   rownames(out) <- NULL
   return(out)
}

# if(F){
#    resamplingGrid(
#       30, 6000, breaks.a=5, breaks.b=5, min.size.diff=15,
#       max.upsample.ratio=5, max.upsample.size=300,
#       min.downsample.ratio=NULL, min.downsample.size=1000,
#       max.imbalance.ratio=15,
#       max.pairs=10, bypass.filter=F
#    )
# }


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
   df, colname.response='response', resample.ratios=NULL, target.sample.sizes=NULL, return.data=T
){
   # if(F){
   #    #df=train[,1:10]
   #    #resample.ratios=c('TRUE'=1.5,'FALSE'=1)
   #    df_raw <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/features/training_labels.txt')
   #    df_raw <- subset(df, in_training_set)
   #
   #    ##
   #    df <- df_raw
   #    rownames(df) <- df$sample
   #    df$response <- df$response=='CNS_Glioma'
   #    colname.response='response'
   #    resample.ratios=c('TRUE'=1.5,'FALSE'=0.5)
   #    target.sample.sizes=NULL
   #
   #    class_weights <- 1 / table(df_raw$response)
   #    class_weights <- class_weights/sum(class_weights)
   #    sample.weights <- class_weights[as.character(df_raw$response)]
   #    names(sample.weights) <- df_raw$sample
   # }

   y <- df[,colname.response]
   indexes <- split(1:nrow(df), y)

   ## Checks --------------------------------
   ##
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

   ##
   if(is.null(names(target.sample.sizes))){
      stop('`target.sample.sizes` must be an integer vector with names (of the sample classes)')
   }
   if(!all(names(target.sample.sizes) %in% names(indexes))){
      stop('Not all class names are present in the names of `target.sample.sizes`')
   }

   ## Main --------------------------------
   indexes_new <- lapply(names(target.sample.sizes), function(i){
      v <- indexes[[i]]
      target_sample_size <- target.sample.sizes[[i]]

      if(length(v)==target_sample_size){ return(v) }

      sample(
         v,
         size=target_sample_size,
         #prob=sample_weights,
         replace=if(target_sample_size > length(v)){ TRUE } else { FALSE }
      )
   })
   #names(indexes_new) <- names(target.sample.sizes)

   indexes_new <- sort(unlist(indexes_new))

   if(return.data){ df[indexes_new,] } else { indexes_new }
}
