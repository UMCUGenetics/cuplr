#' Performs pairwise testing and selects significant features
#'
#' @description For numerical variables, wilcoxon tests are performed. For categorical variables,
#' fisher exact tests are performed. The first factor level is assumed to be the negative outcome,
#' while the other levels are grouped together as the positive outcome. For example,
#' with the factor `as.factor(c('none','loh+pathogenic','deep_deletion'))`, 'none' is considered the
#' negative outcome.
#'
#' When y is a factor (multiclass classification), multiple one-vs-rest pairwise tests (i.e. one for
#' each class label) are performed for each feature. A feature is kept if any of the pairwise tests
#' give a significant pvalue/qvalue.
#'
#' @param x A matrix or dataframe of features
#' @param y A vector of class labels. For binary classification a logical vector. For
#' multiclass classification a factor.
#' @param alternative A vector containing 'two.sided','greater' or 'less', corresponding to each
#' feature
#' @param max.pvalue pvalue threshold for keeping features. Default: 0.01
#' @param min.cliff.delta 0 to 1. Cliff delta threshold for keeping continuous features. Applies to
#' both +ve and -ve values
#' @param min.cramer.v 0 to 1. Cliff delta threshold for keeping categorical features. Applies
#' to both +ve and -ve values
#' @param sel.top.n.features Limit the total number of features that are selected
#' @param min.features Minimum number of features to keep. Prevents outputting no features
#' @param whitelist A character vector of feature names in x to keep, regardless of statistical
#' enrichment
#' @param order.by.pvalue If FALSE, rows of the output dataframe will be ordered by the same
#' order of the features (i.e. colnames) in `x`. If TRUE, it will be sorted by pvalue
#' @param verbose Show progress messages?
#'
#' @return A dataframe containing the pvalue, effect size, and other summary statistics for each
#' feature
#' @export
#'
univarFeatSel <- function(x, y, verbose=F, ...){

   if(is.logical(y)){
      univarFeatSel.logical(x=x, y=y, verbose=verbose, ...)
   } else if(is.factor(y) | is.character(y)) {
      ## Multiclass one vs rest enrichment
      do.call(rbind, lapply(unique(y), function(i){
         #i='Breast'
         if(verbose){ message('## ', i) }
         out <- univarFeatSel(x=x, y=y==i, verbose>=2, ...)
         cbind(class=i, out)
      }))
   } else {
      stop('`y` must be a logical, factor, or character vector')
   }
}

univarFeatSel.logical <- function(
   x, y,
   alternative=NULL,
   max.pvalue=0.01, min.cliff.delta=0.1, min.cramer.v=0.1,
   sel.top.n.features=NULL, min.features=2,
   whitelist=NULL, order.by.pvalue=TRUE,
   verbose=F
){
   if(F){
      base_dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/'
      devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/commonUtils/'))
      devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
      #devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))

      training_data <- read.delim(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.12a_drivers_linx/features/features.txt.gz'))
      training_data <- read.delim(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.09c_originalFeatures/features/features.txt.gz'))

      x=features[,-1]
      x=x[,grep('^rmd',colnames(x),invert=T)]
      y=features[,1]=='Skin'

      max.pvalue=0.01
      min.cliff.delta=0.1
      min.cramer.v=0.1
      min.features=2
      verbose=T

      alternative <- 'two.sided'
      #alternative[ grep('(^purple)|(^rmd)',colnames(x)) ] <- 'two.sided'

      #output.type='new.x'
      sel.top.n.features=NULL
   }

   ## Checks --------------------------------
   if(!is.data.frame(x)){ stop('x must be a dataframe') }
   if(!is.logical(y)){ stop('y must be a logical vector') }
   if(any(sapply(x, is.character))){ stop('characters must be converted to factors') }
   if(is.null(colnames(x))){ stop('x must have colnames') }

   ## Initialize alternative --------------------------------
   if(is.null(alternative)){
      alternative <- rep('two.sided',ncol(x))
   } else if(length(alternative)==1){
      alternative <- rep(alternative,ncol(x))
   } else {
      if( !all(alternative %in% c('two.sided','greater','less')) ){
         stop("`alternative` contains values other than 'two.sided','greater' or 'less'")
      }

      if(length(alternative)!=ncol(x)){
         stop("`alternative` must be the same length as the number of features")
      }
   }

   if(length(names(alternative))==0){
      names(alternative) <- colnames(x)
   }

   ## Convert data types --------------------------------
   if(verbose){ message('Converting factors to logicals (assuming 1st factor level as negative effect)') }
   x <- lapply(x, function(i){
      if(is.factor(i)){ return(as.integer(i)>1) }
      return(i)
   })

   x <- as.data.frame(x, check.names=F)

   ## Wilcox test on numeric features --------------------------------
   is_numeric <- sapply(x, is.numeric)
   x_numeric <- x[,is_numeric, drop=F]
   x_numeric <- as.matrix(x_numeric)

   pvalues_numeric <- numeric()
   cliff_delta <- numeric()
   if(ncol(x_numeric)!=0){
      if(verbose){ message('Performing wilcox tests for numeric features...') }
      pvalues_numeric <- wilcoxTest.data.frame(
         x_numeric[y,,drop=F],
         x_numeric[!y,,drop=F],
         alternative = alternative[colnames(x_numeric)]
      )
      pvalues_numeric[is.na(pvalues_numeric)] <- 1

      if(verbose){ message("Calculating Cliff's delta for numeric features...") }
      cliff_delta <- cliffDelta.data.frame(x_numeric[y,,drop=F], x_numeric[!y,,drop=F])
   }

   ## Fisher test on logical features --------------------------------
   x_logical <- x[,!is_numeric,drop=F]
   x_logical <- as.matrix(x_logical)

   pvalues_logical <- numeric()
   cramer_v <- numeric()
   if(ncol(x_logical)!=0){
      if(verbose){ message('Performing fisher tests for logical features...') }
      x_logical <- as.matrix(as.data.frame(x_logical, check.names=F))
      #table(x_logical[,'fusion.TMPRSS2_ERG'])

      conting <- contingencyMatrix(x_logical, y, use.totals=F)
      pvalues_logical <- fisherTest.data.frame(
         conting,
         alternative = alternative[colnames(x_logical)]
      )

      if(verbose){ message("Calculating Cramer's V for logical features...") }
      cramer_v <- cramerV.data.frame(conting)
   }

   ## Aggregate stats from numeric/logical data --------------------------------
   ## Initialized dataframe
   tests <- data.frame(
      feature=c(colnames(x_numeric), colnames(x_logical)),
      feature_type=c(rep('numeric',ncol(x_numeric)), rep('logical',ncol(x_logical)))
   )
   rownames(tests) <- tests$feature

   ## Add test data
   tests$alternative <- alternative[rownames(tests)]
   tests$pvalue <- c(pvalues_numeric, pvalues_logical)
   tests$cliff_delta <- c(cliff_delta, rep(0,ncol(x_logical)))
   tests$cramer_v <- c(rep(0,ncol(x_numeric)), cramer_v)

   ## Merge cliff delta and cramer's v
   which_not_numeric <- tests$feature_type!='numeric'
   tests$eff_size <- tests$cliff_delta
   tests$eff_size[which_not_numeric] <- tests$cramer_v[which_not_numeric]
   tests$eff_size_metric <- 'cliff_delta'
   tests$eff_size_metric[which_not_numeric] <- 'cramer_v'

   ## Add case/ctrl cohort stats  --------------------------------
   if(verbose){ message("Calculating summary stats...") }
   colMeansTrimmed <- function(m, trim=0.25){
      m <- apply(m,2,sort)

      n <- nrow(m)
      lo <- floor(n * trim) + 1
      hi <- n + 1 - lo

      m <- m[lo:hi,,drop=F]
      colMeans(m)
   }

   calcAvg <- function(x.numeric, x.logical){
      avg_numeric <- numeric()
      if(ncol(x.numeric)!=0){
         avg_numeric <- colMeansTrimmed(x.numeric)
      }

      avg_logical <- numeric()
      if(ncol(x.logical)!=0){
         avg_logical <- colMeans(x.logical)
      }

      c(avg_numeric, avg_logical)
   }

   tests$avg_case <- calcAvg(x_numeric[y,,drop=F], x_logical[y,,drop=F])
   tests$avg_ctrl <- calcAvg(x_numeric[!y,,drop=F], x_logical[!y,,drop=F])
   tests$avg_metric <- 'iqm'
   tests$avg_metric[which_not_numeric] <- 'mean'

   tests <- tests[order(tests$pvalue),]

   ## Post-processing --------------------------------
   ## Select features
   is_pass_feature <- structure(
      rep(FALSE, nrow(tests)),
      names=tests$feature
   )

   is_pass_feature[
      with(tests,{
         pvalue < max.pvalue &
         (
            (alternative %in% c('two.sided','greater') & cliff_delta >=  min.cliff.delta) |
            (alternative %in% c('two.sided','less')    & cliff_delta <= -min.cliff.delta)
         ) |

         (
            (alternative %in% c('two.sided','greater') & cramer_v >=  min.cramer.v) |
            (alternative %in% c('two.sided','less')    & cramer_v <= -min.cramer.v)
         )
      })
   ] <- TRUE

   ## Return feature names or filtered feature matrix
   keep_features <- names(is_pass_feature)[is_pass_feature]
   if(length(keep_features) < min.features){
      keep_features <- names(keep_features)[1:min.features]
   } else if(!is.null(sel.top.n.features)){
      top_n_features <- min(length(keep_features), sel.top.n.features, ncol(x))
      keep_features <- keep_features[ 1:top_n_features ]
   }

   tests$is_pass_feature <- is_pass_feature
   tests$is_keep_feature <- is_pass_feature

   tests$is_keep_feature[ is_pass_feature & !(tests$feature %in% keep_features) ] <- FALSE

   if(!is.null(whitelist)){
      tests$is_keep_feature[ tests$feature %in% whitelist ] <- TRUE
   }

   ## Cleanup output
   NULL -> rownames(tests) -> tests$cliff_delta -> tests$cramer_v

   if(!order.by.pvalue){
      tests <- tests[match(colnames(x), tests$feature),]
   }

   return(tests)
}

