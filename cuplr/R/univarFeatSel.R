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
#' @param v.alternative A vector containing 'two.sided','greater' or 'less', corresponding to each
#' feature
#' @param max.pvalue pvalue threshold for keeping features. Default: 0.01
#' @param min.effect.size 0 to 1 (default: 0.2). Cliff delta threshold for keeping features. Applies
#' to both +ve and -ve cliff delta values
#' @param sel.top.n.features Limit the total number of features that are selected
#' @param output.type Can be 'raw','new.x','features'
#' @param verbose Show progress messages?
#'
#' @return A vector of feature names if return.new.x=TRUE, else a feature matrix with the selected
#' features
#' @export
#'
univarFeatSel <- function(
   x, y,
   v.alternative=NULL,
   max.pvalue=0.01, min.effect.size=0.1, sel.top.n.features=NULL,
   output.type='new.x', verbose=F
){
   if(F){
      base_dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/'
      devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/commonUtils/'))
      devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
      devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))

      training_data <- readFeaturesCuplr(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.09c_originalFeatures/features/features.txt.gz'))
      x=training_data[,-1]
      y=training_data[,1]=='Breast'

      max.pvalue=0.01
      min.effect.size=0.1
      verbose=T

      v.alternative <- rep('greater', ncol(x))
      v.alternative[ grep('(^purple)|(^rmd)',colnames(x)) ] <- 'two.sided'
   }

   ## Checks --------------------------------
   if(!is.data.frame(x)){ stop('x must be a dataframe') }
   if(!is.logical(y)){ stop('y must be a logical vector') }
   if(any(sapply(x, is.character))){ stop('characters must be converted to factors') }
   if(is.null(colnames(x))){ stop('x must have colnames') }
   if(!(output.type %in% c('new.x','raw','features'))){ stop("output.type must be 'new.x', 'raw', or 'features'") }

   ## Initialize alternative --------------------------------
   if(is.null(v.alternative)){
      v.alternative <- rep('two.sided',ncol(x))
   } else {
      if( !all(v.alternative %in% c('two.sided','greater','less')) ){
         stop("`v.alternative` contains values other than 'two.sided','greater' or 'less'")
      }

      if(length(v.alternative)!=ncol(x)){
         stop("`v.alternative` must be the same length as the number of features")
      }
   }
   names(v.alternative) <- colnames(x)

   ## Convert factors to logicals --------------------------------
   if(output.type=='new.x'){ x_raw <- x }

   ## Assume 1st factor level to be the negative effect
   x <- lapply(x, function(i){
      if(is.factor(i)){ return(as.integer(i)>1) }
      return(i)
   })
   x <- as.data.frame(x, check.names=F)

   ## Wilcox test on numeric features --------------------------------
   is_numeric <- sapply(x, is.numeric)
   x_numeric <- x[,is_numeric, drop=F]

   pvalues_numeric <- NULL
   if(ncol(x_numeric)!=0){
      if(verbose){ message('Performing wilcox tests on numerical features...') }
      pvalues_numeric <- wilcoxTest.data.frame(
         x_numeric[y,],
         x_numeric[!y,],
         alternative = v.alternative[colnames(x_numeric)]
      )
      pvalues_numeric[is.na(pvalues_numeric)] <- 1
      names(pvalues_numeric) <- colnames(x_numeric)
   }

   ## Fisher test on logical/categorical features --------------------------------
   x_logical <- x[,!is_numeric,drop=F]

   pvalues_logical <- NULL
   if(ncol(x_logical)!=0){
      if(verbose){ message('Performing fisher tests on logical/categorical features...') }
      x_logical <- as.matrix(as.data.frame(x_logical, check.names=F))
      #table(x_logical[,'fusion.TMPRSS2_ERG'])

      conting <- contingencyMatrix(x_logical, y, use.totals=F)
      pvalues_logical <- fisherTest.data.frame(
         conting,
         alternative = v.alternative[colnames(x_logical)]
      )
      names(pvalues_logical) <- colnames(x_logical)
   }

   pvalues <- c(pvalues_numeric, pvalues_logical)
   pvalues <- pvalues[colnames(x)] ## Preserve original feature order

   ## Post-processing --------------------------------
   tests <- data.frame(
      pvalue=pvalues,
      alternative=v.alternative
   )

   if(verbose){ message('Calculating cliff delta effect size...') }
   tests$cliff_delta <- cliffDelta.data.frame(x[y,], x[!y,])

   tests <- tests[order(tests$pvalue),]

   ## Select features
   tests_ss <- subset(
      tests,
      pvalue < max.pvalue
   )

   if(!is.null(min.effect.size)){
      tests_ss <- subset(
         tests_ss,
         (alternative %in% c('two.sided','greater') & cliff_delta >=  min.effect.size) |
         (alternative %in% c('two.sided','less')    & cliff_delta <= -min.effect.size)
      )
   }

   keep_features <- rownames(tests_ss)
   if(!is.null(sel.top.n.features)){
      keep_features <- keep_features[ 1:min(length(keep_features), sel.top.n.features, ncol(x)) ]
   }

   if(output.type=='new.x'){
      return(x_raw[,keep_features,drop=F])
   }
   if(output.type=='raw'){
      return(tests)
   }

   if(output.type=='features'){
      return(keep_features)
   }
}


