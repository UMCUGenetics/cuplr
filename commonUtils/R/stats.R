#' Perform multiple fisher's exact tests
#'
#' @param m A matrix containing the 4 columns corresponding to a contingency matrix:
#'  xpos: Group 'x' responders
#'  xneg: Group 'x' non responders
#'  ypos: Group 'y' responders
#'  yneg: Group 'y' non responders
#'
#' @param verbose Show progress?
#' @param ... Arguments that can be passed to fisher.test()
#'
#' @return A numeric vector of pvalues
#' @export
#'
fisherTestApply <- function(m, verbose=F, ...){
   #m=df[,3:6]

   if(ncol(m)!=4){ stop('m must have 4 columns: xpos, xneg, ypos, yneg') }
   m <- as.matrix(m)

   if(verbose){
      counter <- 0; pb <- txtProgressBar(max=nrow(m), style=3)
   }

   apply(m, 1, function(i){
      #i=m[1,]
      if(verbose){
         counter <<- counter + 1; setTxtProgressBar(pb, counter)
      }

      m <- matrix(i, nrow=2)
      fisher.test(m, ...)$p.value
   })
}

####################################################################################################
#' Perform multiple wilcox tests
#'
#' @param l A list with pairs of numeric vectors:
#'  list(
#'    first=list(numeric(), numeric()),
#'    second=list(numeric(), numeric()),
#'    third=list(numeric(), numeric())
#'  )
#' @param verbose Show progress?
#' @param ...
#'
#' @return A vector of pvalues
#' @export
#'
wilcoxTestApply <- function(l, verbose=F, ...){
   if(verbose){
      counter <- 0; pb <- txtProgressBar(max=length(l), style=3)
   }
   pvalue <- unlist(lapply(l, function(i){
      if(verbose){ counter <<- counter + 1; setTxtProgressBar(pb, counter) }
      wilcox.test(i[[1]], i[[2]], ...)$p.value
   }))
   pvalue[is.na(pvalue)] <- 1

   return(pvalue)
}


####################################################################################################
#' Calculate effect size
#'
#' @param xpos Group 'x' responders
#' @param xneg Group 'x' non responders
#' @param ypos Group 'y' responders
#' @param yneg Group 'y' non responders
#' @param m A matrix containing the xpos, xneg, ypos, yneg (in this order)
#'
#' @description Effect size calculation is based on odds ratio = (xpos/xneg) / (ypos/yneg). Odds
#' ratios < 1 are scaled to become negative values.
#'
#' @return A numeric vector
#' @export
#'
effectSize <- function(xpos=NULL, xneg=NULL, ypos=NULL, yneg=NULL, m=NULL){

   if(is.matrix(xpos) | is.data.frame(xpos)){
      stop('`xpos` is a matrix/dataframe. Use the `m` argument instead')
   }
   if(!is.null(m)){
      xpos <- m[,1]
      xneg <- m[,2]
      ypos <- m[,3]
      yneg <- m[,4]
   }

   ## Calc preliminary odds ratios
   x <- xpos/xneg
   y <- ypos/yneg

   xy <- x/y
   yx <- y/x

   out <- xy

   ## Scale odds ratios < 1 to be negative
   out[out<1] <- NA
   out[is.na(out)] <- -(yx[is.na(out)])

   ## Deal with 0/0 cases
   out[is.na(out)] <- 0

   return(out)
}

continuousToContingency <- function(x, y, mid.func=median){
   #x=test_pairs[[1]][[1]]
   #y=test_pairs[[1]][[2]]

   xmid <- mid.func(x)
   ymid <- mid.func(y)

   c(
      xpos=sum(x >  ymid),
      xneg=sum(x <= ymid),
      ypos=sum(y >  xmid),
      yneg=sum(y <= xmid)
   )
}






