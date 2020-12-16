#' Calculate Cliff's delta values
#'
#' @rdname cliffDelta
#'
#' @description Calculate effect sizes using Cliff's delta values for pairwise continuous variable
#' comparisons (e.g. comparisons that would be done by a wilcox test or t-test)
#'
#' @param x A numeric vector for the 1st group of observations. Alternatively, a matrix, where
#' comparisons will be performed by row.
#' @param y Same as x but for the 2nd group of observations.
#'
#' @return A numeric vector of Cliff's delta values
#' @export
#'
cliffDelta <- function (x, ...) {
   UseMethod("cliffDelta", x)
}

Rcpp::cppFunction('double cliffd(NumericVector x, NumericVector y){
   int len_x = x.size();
   int len_y = y.size();

   int sign_sum = 0;
   for(int i = 0; i < len_x; i++){
      for(int j = 0; j < len_y; j++){
         if(x[i] < y[j]){
            sign_sum--;
         } else if(x[i] > y[j]){
            sign_sum++;
         }
      }
   }

   return sign_sum / ((double)len_x*len_y);
}
')

# set.seed(1)
# x <- sample(1:10,6,replace=T)
# y <- sample(1:10,5,replace=T)
#
# x <- m_true[,1]
# y <- m_false[,1]
#
# ## R version
# signs <- sign(outer(x, y, FUN="-")); sum(signs, na.rm=T) / length(signs)
#
# ## Cpp version
# cliffd(x,y)

#' @rdname cliffDelta
#' @method cliffDelta default
cliffDelta.default <- function(x,y){
   #x <- m_true[,1]
   #y <- m_false[,1]
   if(!is.numeric(x) | !is.numeric(y)){ stop('x and y must be numeric matrices') }
   cliffd(x,y)
}

#' @rdname cliffDelta
#' @method cliffDelta matrix
cliffDelta.matrix  <- function(x, y){
   #x=m_true
   #y=m_false

   if(ncol(x)!=ncol(y)){ stop('x and y must have the sample number of columns') }
   x <- as.matrix(x); dimnames(x) <- NULL
   y <- as.matrix(y); dimnames(y) <- NULL

   if(!is.numeric(x) | !is.numeric(y)){ stop('x and y must be numeric matrices') }

   # ## R implementation
   # outer_sign_sum <- unlist(lapply(1L:ncol(x), function(i){
   #    #i=1
   #    sum(
   #       outer(x[,i], y[,i], FUN=function(a,b){ sign(a-b) }),
   #       na.rm=T
   #    )
   # }), use.names=F)
   # outer_sign_sum / (nrow(x)*nrow(y))

   ## Cpp implementation
   unlist(lapply(1L:ncol(x), function(i){
      #i=1
      cliffd(x[,i],y[,i])
   }), use.names=F)
}

#' @rdname cliffDelta
#' @method cliffDelta dataframe
cliffDelta.data.frame <- cliffDelta.matrix

####################################################################################################
#' Calculate odds ratio
#'
#' @param xpos Group 'x' responders
#' @param xneg Group 'x' non responders
#' @param ypos Group 'y' responders
#' @param yneg Group 'y' non responders
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
oddsRatio.default <- function(xpos,xneg,ypos,yneg, trans.method='none', logistic.growth=0.2){
   # if(F){
   #    xpos=c(50, 0,50, 10,1000)
   #    xneg=c(60, 0,0,  20,1001)
   #    ypos=c(30, 0,50, 20,2)
   #    yneg=c(200,0,100,0, 500)
   # }

   ## Calc preliminary odds ratios
   x <- xpos/xneg
   y <- ypos/yneg

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
oddsRatio.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   oddsRatio.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname oddsRatio
#' @method oddsRatio dataframe
oddsRatio.data.frame <- oddsRatio.matrix



