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
#' @param verbose Show progress messages?
#'
#' @return A vector of feature names if return.new.x=TRUE, else a feature matrix with the selected
#' features
#' @export
#'
univarFeatSel <- function(
   x, y,
   alternative=NULL,
   max.pvalue=0.01, min.cliff.delta=0.1, min.cramer.v=0.1,
   sel.top.n.features=NULL, min.features=2,
   whitelist=NULL,
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
   #if(!(output.type %in% c('new.x','raw','features'))){ stop("output.type must be 'new.x', 'raw', or 'features'") }

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

   ## Convert factors to logicals --------------------------------
   #if(output.type=='new.x'){ x_raw <- x }

   if(verbose){ message('Converting factors to logicals (assuming 1st factor level as negative effect)') }
   x <- lapply(x, function(i){
      if(is.factor(i)){ return(as.integer(i)>1) }
      return(i)
   })
   x <- as.data.frame(x, check.names=F)

   ## Wilcox test on numeric features --------------------------------
   is_numeric <- sapply(x, is.numeric)
   x_numeric <- x[,is_numeric, drop=F]

   pvalues_numeric <- numeric()
   cliff_delta <- numeric()
   if(ncol(x_numeric)!=0){
      if(verbose){ message('Performing wilcox tests for numeric features...') }
      pvalues_numeric <- wilcoxTest.data.frame(
         x_numeric[y,],
         x_numeric[!y,],
         alternative = alternative[colnames(x_numeric)]
      )
      pvalues_numeric[is.na(pvalues_numeric)] <- 1

      if(verbose){ message("Calculating Cliff's delta for numeric features...") }
      cliff_delta <- cliffDelta.data.frame(x_numeric[y,], x_numeric[!y,])
   }

   ## Fisher test on logical features --------------------------------
   x_logical <- x[,!is_numeric,drop=F]

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
   tests <- data.frame(
      feature_type=c(rep('numeric',ncol(x_numeric)), rep('logical',ncol(x_logical))),
      pvalue=c(pvalues_numeric, pvalues_logical),
      cliff_delta=c(cliff_delta, rep(0,ncol(x_logical))),
      cramer_v=c(rep(0,ncol(x_numeric)), cramer_v),
      row.names=c(colnames(x_numeric), colnames(x_logical))
   )

   tests <- tests[colnames(x),] ## Preserve original feature order
   tests$alternative <- alternative[rownames(tests)]

   tests <- data.frame(feature=rownames(tests), tests, row.names=NULL)
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

   return(tests)
}


