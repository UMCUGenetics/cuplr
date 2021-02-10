#' Calculate Cliff's delta values (for continuous data)
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
#' @export
cliffDelta.default <- function(x,y){
   #x <- m_true[,1]
   #y <- m_false[,1]
   if(!is.numeric(x) | !is.numeric(y)){ stop('x and y must be numeric matrices') }
   cliffd(x,y)
}

#' @rdname cliffDelta
#' @method cliffDelta matrix
#' @export
cliffDelta.matrix  <- function(x, y){
   #x=m_true
   #y=m_false

   if(ncol(x)!=ncol(y)){ stop('x and y must have the sample number of columns') }
   x <- as.matrix(x); dimnames(x) <- NULL
   y <- as.matrix(y); dimnames(y) <- NULL

   if(!is.numeric(x) & !is.logical(x)){ stop('x must be a numeric or logical matrix') }
   if(!is.numeric(y) & !is.logical(y)){ stop('y must be a numeric or logical matrix') }

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
#' @method cliffDelta data.frame
#' @export
cliffDelta.data.frame <- cliffDelta.matrix


####################################################################################################
#' Calculate odds ratio
#'
#' @rdname cramerV
#'
#' @description Calculate effect sizes using Cramer's V values for pairwise categorical variable
#' comparisons (e.g. comparisons that would be done by a fisher's exact test)
#'
#' @param case.true Case group responders
#' @param case.false Case group non responders
#' @param ctrl.true Control group responders
#' @param ctrl.false Control group non responders
#' @param show.sign Add `-` to values where effect decreases?
#'
#' @return A numeric vector
#' @export
#'
cramerV <- function (x, ...) {
   UseMethod("cramerV", x)
}

#' @rdname cramerV
#' @method cramerV default
#' @export
cramerV.default <- function(case.true, case.false, ctrl.true, ctrl.false, show.sign=T){
   # if(F){
   #    case.true= c(290 ,60)
   #    case.false=c(289 ,247)
   #    ctrl.true= c(3   ,0)
   #    ctrl.false=c(5034,5309)
   # }

   ## The 2x2 contigency matrix looks like this
   ## case.true   ctrl.true
   ## case.false  ctrl.false

   ## Calculate expect values --------------------------------
   ## Based on this tutorial: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
   ## From chisq.test()
   # x <- matrix(c(case.true[1], case.false[1], ctrl.true[1], ctrl.false[1]), nrow=2)
   # n <- sum(x)
   # sr <- rowSums(x)
   # sc <- colSums(x)
   # E <- outer(sr, sc, "*")/n

   ## Vectorized version
   observed <- cbind(case.true, case.false, ctrl.true, ctrl.false)
   n <- rowSums(observed)

   rowsums1 <- case.true  + ctrl.true
   rowsums2 <- case.false + ctrl.false
   colsums1 <- case.true  + case.false
   colsums2 <- ctrl.true  + ctrl.false

   expected <- cbind(
      case.true  = rowsums1*colsums1,
      case.false = rowsums2*colsums1,
      ctrl.true  = rowsums1*colsums2,
      ctrl.false = rowsums2*colsums2
   ) / n

   ## Chi-squared statistic
   chi2 <- rowSums(
      (observed - expected)^2 / expected
   )

   ## Calculate Cramers V --------------------------------
   ## Based on this tutorial: https://www.real-statistics.com/chi-square-and-f-distributions/effect-size-chi-square/
   ## V = sqrt(chi2/(n*df))
   ## where df* = min(r – 1, c – 1) and r = the number of rows and c = the number of columns in the contingency table.
   V <- sqrt(chi2/n)
   V[is.na(V)] <- 0

   if(show.sign){
      ## Add sign
      observed_1 <- observed + 1
      ratio_case <- observed_1[,'case.true'] / observed_1[,'case.false']
      ratio_ctrl <- observed_1[,'ctrl.true'] / observed_1[,'ctrl.false']
      log_odds <- log2(ratio_case / ratio_ctrl)
      V[log_odds<0] <- -V[log_odds<0]
   }

   return(V)
}

#' @rdname cramerV
#' @method cramerV matrix
#' @export
cramerV.matrix <- function(m){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   cramerV.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname cramerV
#' @method cramerV data.frame
#' @export
cramerV.data.frame <- cramerV.matrix


####################################################################################################
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



