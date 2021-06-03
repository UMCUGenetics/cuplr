#' Calculate odds ratio
#'
#' @param case.true Case group responders
#' @param case.false Case group non responders
#' @param ctrl.true Control group responders
#' @param ctrl.false Control group non responders
#' @param trans.method Can be 'none', 'linear' or 'logistic'
#' @param logistic.growth Steepness/gradient of the logistic function when trans.method=='logistic'
#'
#' @description Effect size calculation is based on odds ratio = (xpos/xneg) / (ypos/yneg). When
#' trans.method=='linear' odds ratios < 1 are scaled to become negative values. When
#' trans.method=='logistic' the linearly scaled odds ratios are logistically transformed with the
#' equation 2 / (1 + exp(1)^(-logistic.growth*OR) ) - 1
#'
#' @return A numeric vector
#' @export
#'
oddsRatio <- function (x, ...) {
   UseMethod("oddsRatio", x)
}

#' @rdname oddsRatio
#' @method oddsRatio default
oddsRatio.default <- function(case.true, case.false, ctrl.true, ctrl.false, trans.method='none', logistic.growth=0.2){
   # if(F){
   #    xpos=c(50, 0,50, 10,1000)
   #    xneg=c(60, 0,0,  20,1001)
   #    ypos=c(30, 0,50, 20,2)
   #    yneg=c(200,0,100,0, 500)
   # }

   ## Calc preliminary odds ratios
   x <- case.true/case.false
   y <- ctrl.true/ctrl.false

   or <- x/y
   or_inv <- y/x

   # if(trans.method=='chinn'){
   #    ## Chinn S., 2000
   #    ## https://www.real-statistics.com/chi-square-and-f-distributions/effect-size-chi-square/
   #    out <- log(or)/1.81379936423422 ## 1.81 = pi/sqrt(3)
   #    out[is.na(out)] <- 0
   #    return(out)
   # }

   ## Scale odds ratios < 1 to be negative
   out <- or
   out[out<1] <- NA
   out[is.na(out)] <- -(or_inv[is.na(out)])

   ## Deal with 0/0 cases
   out[is.na(out)] <- 0

   if(trans.method=='linear'){
      return(out)
   }

   if(trans.method=='logistic'){
      # df <- data.frame(x=seq(-10,10,by=0.5))
      # k=0.2
      # df$y <-  2 / (1 + exp(1)^(-k*df$x) ) - 1
      # df
      return(
         2 / (1 + exp(1)^(-logistic.growth*out) ) - 1
      )
   }

   return(or)
}

#' @rdname oddsRatio
#' @method oddsRatio matrix
#' @export
oddsRatio.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   oddsRatio.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname oddsRatio
#' @method oddsRatio data.frame
#' @export
oddsRatio.data.frame <- oddsRatio.matrix



