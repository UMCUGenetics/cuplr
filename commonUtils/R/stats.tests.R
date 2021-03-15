#' Create multiple contingency matrices
#'
#' @param x A logical matrix
#' @param y A logical vector indicating the case/control group (i.e. response)
#' @param use.totals Should case/control totals be used instead of case.false and ctrl.false?
#' @param show.warnings Show warnings?
#'
#' @return A matrix with the columns: case.true, case.false, ctrl.true, ctrl.false. Each row is
#' thus one contingency matrix
#' @export
#'
contingencyMatrix <- function(x, y, use.totals=FALSE, show.warnings=TRUE){

   if(!is.logical(x)){ stop('x must be a logical matrix or vector') }
   if(!is.logical(y)){ stop('y must be a logical vector') }

   if(!is.matrix(x)){ x <- matrix(x, ncol=1) }

   case <- x[y,,drop=F]
   ctrl <- x[!y,,drop=F]

   if(anyNA(x)){
      if(show.warnings){
         warning('`x` contains NAs. These were excluded from the contingency counts')
      }
      case.NAs <- apply(case, 2, function(i){ sum(is.na(i)) })
      ctrl.NAs <- apply(ctrl, 2, function(i){ sum(is.na(i)) })

      case[is.na(case)] <- FALSE
      ctrl[is.na(ctrl)] <- FALSE
   } else {
      case.NAs <- 0
      ctrl.NAs <- 0
   }

   case.true <- colSums(case)
   case.false <- nrow(case) - case.true - case.NAs

   ctrl.true <- colSums(ctrl)
   ctrl.false <- nrow(ctrl) - ctrl.true - ctrl.NAs

   if(!use.totals){
      cbind(
         case.true, case.false,
         ctrl.true, ctrl.false
      )
   } else {
      cbind(
         case.true,
         case.total=nrow(case) - case.NAs, ## Use case totals instead
         ctrl.true,
         case.total=nrow(ctrl) - ctrl.NAs ## Use ctrl totals instead
      )
   }
}

#' Fast fisher's exact tests
#'
#' @rdname fisherTest
#'
#' @param case.true Case group (treated) responders
#' @param case.false Case group (treated) non responders
#' @param ctrl.true Control group (untreated) responders
#' @param ctrl.false Control group (untreated) non responders
#' @param m A matrix containing the 4 columns corresponding to a contingency matrix:
#'
#' @param verbose Show progress?
#' @param ... Arguments that can be passed to fisherTest.default()
#'
#' @return A numeric vector of pvalues
#' @export
#'
fisherTest <- function (x, ...) {
   UseMethod("fisherTest", x)
}

#' @rdname fisherTest
#' @method fisherTest default
#' @export
fisherTest.default <- function(
   case.true, case.false, ctrl.true, ctrl.false,
   alternative='two.sided', verbose=F
){

   ## Based on:
   ## https://stats.stackexchange.com/questions/454248/fisher-test-alternative-problem-in-r
   x <- case.true
   m <- case.true + case.false
   n <- ctrl.true + ctrl.false
   k <- case.true + ctrl.true

   if(verbose){
      pb <- txtProgressBar(max=length(case.true), style=3L)
      counter <- 0L
   }
   unlist(
      Map(function(x,m,n,k,alternative){

         if(verbose){
            counter <<- counter + 1L
            setTxtProgressBar(pb, counter)
         }

         a <- 0L:min(m, k)
         prob <- dhyper(a, m, n, k)

         sel_probs <- switch(
            alternative,
            #two.sided=prob[prob<=prob[a==x]],
            two.sided=prob[prob<=prob[x+1L]],
            greater=prob[a>=x],
            less=prob[a<=x]
         )

         sum(sel_probs)
      }, x,m,n,k,alternative),

      use.names=F
   )
}

#' @rdname fisherTest
#' @method fisherTest matrix
#' @export
fisherTest.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   fisherTest.default(
      m[,1L], m[,2L], m[,3L], m[,4L],
      ...
   )
}

#' @rdname fisherTest
#' @method fisherTest data.frame
#' @export
fisherTest.data.frame <- fisherTest.matrix


#' Fast chi-squared test for 2x2 matrices (similar to a fisher test)
#'
#' @rdname fisherChi2
#'
#' @param case.true Group 'x' responders
#' @param case.false Group 'x' non responders
#' @param ctrl.true Group 'y' responders
#' @param ctrl.false Group 'y' non responders
#' @param m A matrix containing the 4 columns corresponding to a contingency matrix
#'
#' @param verbose Show progress?
#' @param ... Arguments that can be passed to fisherChi2.default()
#'
#' @return A numeric vector of pvalues or chi-squared statistics
#' @export
#'
fisherChi2 <- function (x, ...) {
   UseMethod("fisherChi2", x)
}

#' @rdname fisherChi2
#' @method fisherChi2 default
#' @export
fisherChi2.default <- function(case.true, case.false, ctrl.true, ctrl.false, correct=TRUE, return.statistic=FALSE){
   ## Adapted from stats::chisq.test()
   # if(F){
   #    #x=matrix(c(50,60,30,200),nrow=2)
   #    case.true=c(50,50,10)
   #    case.false=c(60,60,20)
   #    ctrl.true=c(30,30,200)
   #    ctrl.false=c(200,200,1000)
   # }

   m <- cbind(case.true, case.false, ctrl.true, ctrl.false)
   n <- rowSums(m)

   sums <- data.frame(
      row1=case.true + ctrl.true,
      row2=case.false + ctrl.false,
      col1=case.true + case.false,
      col2=ctrl.true + ctrl.false
   )

   E <- do.call(rbind, Map(function(row1,row2,col1,col2, n){
      as.vector(outer(c(row1,row2), c(col1,col2), "*")/n)
   }, sums$row1,sums$row2,sums$col1,sums$col2, n))

   m_minus_E <- abs(m - E)

   YATES <- 0
   if(correct){
      YATES <- apply( cbind(0.5, m_minus_E),1,min )
   }

   STATISTIC <- rowSums((m_minus_E - YATES)^2/E)
   if(return.statistic){ return(STATISTIC) }

   pchisq(STATISTIC, 1, lower.tail = FALSE)
}

#' @rdname fisherChi2
#' @method fisherChi2 matrix
#' @export
fisherChi2.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   fisherChi2.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname fisherChi2
#' @method fisherChi2 data.frame
#' @export
fisherChi2.data.frame <- fisherChi2.matrix


####################################################################################################
## Adapted from https://github.com/karoliskoncevicius/matrixTests/

## Helper functions ========================
if(T){
   ## Asserts
   assert_numeric_mat_or_vec <- function(x) {
      name <- as.character(substitute(x))
      if(is.null(x) || !is.numeric(x) | !(is.matrix(x) | is.vector(x)))
         stop(paste0('"', name, '"', ' must be a numeric matrix or vector'))
   }

   assert_logical_vec_length <- function(x, ...) {
      name   <- as.character(substitute(x))
      lens   <- unlist(list(...))
      lnames <- as.character(substitute(list(...)))[-1]
      lnames <- paste(lnames, collapse=' or ')
      if(!(length(x) %in% lens) | !is.logical(x) | (NCOL(x) > 1 & NROW(x) > 1))
         stop(paste0('"', name, '"', ' must be a logical vector with length ', lnames))
   }

   assert_character_vec_length <- function(x, ...) {
      name   <- as.character(substitute(x))
      lens   <- unlist(list(...))
      lnames <- as.character(substitute(list(...)))[-1]
      lnames <- paste(lnames, collapse=' or ')
      if(!(length(x) %in% lens) | !is.character(x) | (NCOL(x) > 1 & NROW(x) > 1))
         stop(paste0('"', name, '"', ' must be a character vector with length ', lnames))
   }

   assert_numeric_vec_length <- function(x, ...) {
      name   <- as.character(substitute(x))
      lens   <- unlist(list(...))
      lnames <- as.character(substitute(list(...)))[-1]
      lnames <- paste(lnames, collapse=' or ')
      if(!(length(x) %in% lens) | !is.numeric(x) | (NCOL(x) > 1 & NROW(x) > 1))
         stop(paste0('"', name, '"', ' must be a numeric vector with length ', lnames))
   }

   assert_all_in_set <- function(x, vals) {
      name <- as.character(substitute(x))
      vnames <- paste(vals, collapse=", ")
      if(is.null(x) | !all(x %in% vals))
         stop(paste0('all "', name, '" values must be in: ', vnames))
   }

   assert_all_in_open_interval <- function(x, min, max) {
      name <- as.character(substitute(x))
      if(is.null(x) | any(anyNA(x) | x<=min | x>=max))
         stop(paste0('all "', name, '" values must be greater than ', min, ' and lower than ', max))
   }

   assert_equal_nrow <- function(x, y) {
      namex <- as.character(substitute(x))
      namey <- as.character(substitute(y))
      if(nrow(x) != nrow(y))
         stop(paste0('"', namex, '" and "', namey, '" must have the same number of rows'))
   }

   ## General
   showWarning <- function(isWarning, err) {
      if(any(isWarning, na.rm=TRUE)) {
         parentFun <- deparse(as.list(sys.call(-1))[[1]])
         grandFun  <- as.list(sys.call(-2))
         if(length(grandFun) > 0) {
            grandFun <- deparse(grandFun[[1]])
            if(grandFun %in% getNamespaceExports("matrixTests")) {
               parentFun <- grandFun
            }
         }
         pref <- "row"
         if(grepl("^col_", parentFun)) pref <- "column"
         n <- sum(isWarning, na.rm=TRUE)
         i <- match(TRUE, isWarning)
         err <- paste0(parentFun, ": ", n, ' of the ', pref, 's ', err, ".",
                       '\nFirst occurrence at ', pref, ' ', i
         )
         warning(err, call.=FALSE)
      }
   }

   rowTies <- function(x) {
      dupRows <- apply(x, 1, anyDuplicated, incomparables=NA) != 0
      if(any(dupRows)) {
         dups <- matrix(FALSE, nrow=nrow(x), ncol=ncol(x))
         dups[dupRows,] <- t(apply(x[dupRows,,drop=FALSE], 1, duplicated, incomparables=NA))
         dups <- cbind(which(dups, arr.ind=TRUE), val=x[dups])
         dups <- dups[order(dups[,1], dups[,3]),,drop=FALSE]

         sp <- split(dups[,3], dups[,1])
         sp <- lapply(sp, function(x) rle(x)$length+1)
         cl <- lapply(sp, seq_along)

         dups <- cbind(rep(unique(dups[,1]), lengths(sp)),
                       unlist(cl),
                       unlist(sp)
         )

         res <- matrix(0L, nrow=nrow(x), ncol=max(dups[,2]))
         res[dups[,1:2,drop=FALSE]] <- dups[,3]
      } else {
         res <- matrix(0L, nrow=nrow(x), ncol=1)
      }
      res
   }

   ## Tests
   do_wilcox_2_exact <- function(stat, nx, ny, alt) {
      res <- rep(NA_integer_, length(stat))

      case <- stat > (nx*ny*0.5)


      inds <- alt=="two.sided" & case
      if(any(inds)) {
         res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
         res[inds] <- pmin(2*res[inds], 1)
      }

      inds <- alt=="two.sided" & !case
      if(any(inds)) {
         res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
         res[inds] <- pmin(2*res[inds], 1)
      }

      inds <- alt=="greater"
      if(any(inds)) {
         res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
      }

      inds <- alt=="less"
      if(any(inds)) {
         res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
      }

      res
   }

   do_wilcox_2_approx <- function(stat, nx, ny, alt, nties, correct) {
      res <- rep(NA_integer_, length(stat))

      z <- stat - nx*ny*0.5
      correction <- rep(0, length(stat))
      correction[correct & alt=="two.sided"] <- sign(z[correct & alt=="two.sided"]) * 0.5
      correction[correct & alt=="greater"]   <- 0.5
      correction[correct & alt=="less"   ]   <- -0.5
      z <- z - correction

      sigma <- sqrt((nx*ny/12) * ((nx+ny+1) - rowSums(nties^3 - nties, na.rm=TRUE) / ((nx+ny) * (nx+ny-1))))
      z <- z/sigma


      inds <- alt=="two.sided"
      if(any(inds)) {
         res[inds] <- 2 * pmin(stats::pnorm(z[inds]), stats::pnorm(z[inds], lower.tail=FALSE))
      }

      inds <- alt=="greater"
      if(any(inds)) {
         res[inds] <- stats::pnorm(z[inds], lower.tail=FALSE)
      }

      inds <- alt=="less"
      if(any(inds)) {
         res[inds] <- stats::pnorm(z[inds])
      }

      res
   }

}

## Main ========================
#' Perform multiple wilcox tests
#'
#' @rdname wilcoxTest
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param mu true values of the location shift for the null hypothesis.
#' A single number or numeric vector with values for each observation.
#' @param exact logical or NA (default) indicator whether an exact p-value
#' should be computed (see Details).
#' A single value or a logical vector with values for each observation.
#' @param correct logical indicator whether continuity correction should be
#' applied in the cases where p-values are obtained using normal approximation.
#' A single value or logical vector with values for each observation.
#' @param pvalue.only Only return pvalues?
#'
#' @return If pvalue.only==TRUE, a numeric vector of pvalues
#' Else, a data.frame where each row contains the results of a wilcoxon test
#' performed on the corresponding row/column of x.
#' The columns will vary depending on the type of test performed.\cr\cr
#' They will contain a subset of the following information:\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 4. obs.paired - number of paired observations (present in x and y)\cr
#' 5. statistic - Wilcoxon test statistic\cr
#' 6. pvalue - p-value\cr
#' 7. alternative - chosen alternative hypothesis\cr
#' 8. location.null - location shift of the null hypothesis\cr
#' 9. exact - indicates if exact p-value was computed\cr
#' 10. correct - indicates if continuity correction was performed
#'
wilcoxTest <- function (x, ...) {
   UseMethod("wilcoxTest", x)
}

#' @rdname wilcoxTest
#' @method wilcoxTest matrix
#' @export
wilcoxTest.matrix <- function(
   x, y, alternative="two.sided", mu=0, exact=NA, correct=TRUE,
   pvalue.only=TRUE
){
   ## Function from row_wilcoxon_twosample()
   ## Translate to col_wilcoxon_twosample()
   ## Therefore need to transpose matrices
   x <- t(x)
   y <- t(y)

   force(x)
   force(y)

   # if(is.vector(x))
   #    x <- matrix(x, nrow=1)
   # if(is.vector(y))
   #    y <- matrix(y, nrow=1)

   # if (is.data.frame(x) && all(sapply(x, is.numeric)))
   #    x <- data.matrix(x)
   # if(is.data.frame(y) && all(sapply(y, is.numeric)))
   #    y <- data.matrix(y)

   assert_numeric_mat_or_vec(x)
   assert_numeric_mat_or_vec(y)

   if(nrow(y)==1L & nrow(x)>1L) {
      y <- matrix(y, nrow=nrow(x), ncol=ncol(y), byrow=TRUE)
   }

   assert_equal_nrow(x, y)

   if(length(alternative)==1)
      alternative <- rep(alternative, length.out=nrow(x))
   assert_character_vec_length(alternative, 1, nrow(x))

   choices <- c("two.sided", "less", "greater")
   alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
   assert_all_in_set(alternative, choices)

   if(length(mu)==1)
      mu <- rep(mu, length.out=nrow(x))
   assert_numeric_vec_length(mu, 1, nrow(x))
   assert_all_in_open_interval(mu, -Inf, Inf)

   if(length(exact)==1)
      exact <- rep(exact, length.out=nrow(x))
   assert_logical_vec_length(exact, 1, nrow(x))
   assert_all_in_set(exact, c(TRUE, FALSE, NA))

   if(length(correct)==1)
      correct <- rep(correct, length.out=nrow(x))
   assert_logical_vec_length(correct, 1, nrow(x))
   assert_all_in_set(correct, c(TRUE, FALSE))


   hasinfx <- is.infinite(x)
   x[hasinfx] <- NA
   hasinfx <- rowSums(hasinfx) > 0

   hasinfy <- is.infinite(y)
   y[hasinfy] <- NA
   hasinfy <- rowSums(hasinfy) > 0

   nxs  <- rep.int(ncol(x), nrow(x)) - matrixStats::rowCounts(is.na(x))
   nys  <- rep.int(ncol(y), nrow(y)) - matrixStats::rowCounts(is.na(y))

   naexact <- is.na(exact)
   exact[naexact] <- (nxs[naexact] < 50) & (nys[naexact] < 50)

   r <- matrixStats::rowRanks(cbind(x - mu, y), ties.method="average")

   statistic <- rowSums(r[,seq_len(ncol(x)),drop=FALSE], na.rm=TRUE) - nxs * (nxs + 1)*0.5

   nties   <- rowTies(r)
   hasties <- rowSums(nties>0) > 0

   wres <- rep(NA_integer_, nrow(x))
   inds <- exact & !hasties
   wres[inds]  <- do_wilcox_2_exact(statistic[inds], nxs[inds], nys[inds], alternative[inds])
   wres[!inds] <- do_wilcox_2_approx(
      statistic[!inds], nxs[!inds], nys[!inds], alternative[!inds],
      nties[!inds,,drop=FALSE], correct[!inds]
   )

   w1 <- hasinfx
   showWarning(w1, 'had infinite "x" observations that were removed')

   w2 <- hasinfy
   showWarning(w2, 'had infinite "y" observations that were removed')

   w3 <- nxs < 1
   showWarning(w3, 'had less than 1 remaining finite "x" observation')

   w4 <- nys < 1
   showWarning(w4, 'had less than 1 remaining finite "y" observation')

   w5 <- exact & hasties
   showWarning(w5, 'had ties: cannot compute exact p-values with ties')

   statistic[w3 | w4] <- NA

   exact <- exact & !hasties
   correct <- correct & !exact


   rnames <- rownames(x)
   if(!is.null(rnames)) rnames <- make.unique(rnames)

   if(pvalue.only){ return(wres) }

   data.frame(
      obs.x=nxs, obs.y=nys, obs.tot=nxs+nys, statistic=statistic,
      pvalue=wres, alternative=alternative, location.null=mu,
      exact=exact, corrected=correct,
      stringsAsFactors=FALSE, row.names=rnames
   )
}

#' @rdname wilcoxTest
#' @method wilcoxTest default
#' @export
wilcoxTest.default <- function(x, y, ...){
   x <- matrix(x, ncol=1L)
   y <- matrix(y, ncol=1L)
   wilcoxTest.matrix(x,y, ...)
}

#' @rdname wilcoxTest
#' @method wilcoxTest data.frame
#' @export
wilcoxTest.data.frame <- wilcoxTest.matrix

#' @rdname wilcoxTest
#' @method wilcoxTest default
#' @export
wilcoxTest.default <- function(x, y, ...){
   x <- matrix(x, ncol=1L)
   y <- matrix(y, ncol=1L)
   wilcoxTest.matrix(x,y, ...)
}



