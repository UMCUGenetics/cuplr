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
#' @param max.qvalue qvalue threshold for keeping features
#' @param max.pvalue pvalue threshold for keeping features. Overrides max.qvalue
#' @param sel.top.n.features Limit the total number of features that are selected
#' @param return.new.x If TRUE, returns a feature matrix with the selected features. Else, a
#' character vector of the selected features
#' @param verbose Show progress messages?
#'
#' @return A vector of feature names if return.new.x=TRUE, else a feature matrix with the selected
#' features
#' @export
#'
univarFeatSel <- function(
   x, y,
   max.qvalue=0.01, max.pvalue=NULL, sel.top.n.features=NULL,
   return.new.x=T, verbose=F
){
   # colname.response='response'
   # x <- df[,colnames(df)!=colname.response]
   # y <- as.factor(df[,colname.response])
   # y <- y=='Prostate'

   if( !(is.logical(y) | is.factor(y)) ){
      stop('`y` must be a logical or factor')
   }
   #if(is.logical(y)){ y <- factor(y,c('TRUE','FALSE')) }

   main <- function(v, y.logical){
      #y.logical=y

      ## Numeric data: wilcox text
      if(is.numeric(v)){
         #v=x$viral_ins.Hepatitis_C_virus
         #v=x$rmd.14q_107
         v_split <- split(v, y.logical)
         wilcox.test(v_split[['TRUE']], v_split[['FALSE']])$p.value

      ## Categorical data: fisher test
      } else {
         #v=ifelse(x$purple.gender,'male','female')
         #v=unname(m[,'AR'])
         #v=x$gene_def.VHL
         #y=metadata[ match(rownames(m), metadata$sample),'cancer_type' ]
         #y.logical <- y=='Prostate'

         neg_category <- levels(as.factor(v))[1]
         #v=rep('0;none',length(v))

         fisher.test(
            matrix(
               c(
                  sum(v!=neg_category & y.logical), sum(y.logical),
                  sum(v!=neg_category & !y.logical), length(y.logical)
               ),
               nrow=2
            )
         )$p.value
      }
   }

   ## Binary classification
   if(is.logical(y)){
      if(verbose){ counter <- 0 }
      p_values <- unlist(lapply(as.data.frame(x), function(i){
         if(verbose){
            counter <<- counter + 1
            message('[',counter,'] ', colnames(x)[[counter]] )
         }
         main(i, y)
      }))
      p_values <- sort(p_values)
      q_values <- p.adjust(p_values, method='bonferroni')

      if(!is.null(max.pvalue)){
         keep_features <- names(p_values)[ p_values < max.pvalue ]
      } else {
         keep_features <- names(q_values)[ q_values < max.qvalue ]
      }

      if(!is.null(sel.top.n.features)){ keep_features <- keep_features[1:sel.top.n.features] }
      keep_features <- na.exclude(keep_features)

   ## Multiclass classification
   } else {
      y_logicals <- lapply(levels(y), function(i){ y==i })
      names(y_logicals) <- levels(y)

      if(verbose){ counter <- 0 }
      m_p_values <- do.call(cbind, lapply(y_logicals, function(y_logical){
         if(verbose){
            counter <<- counter + 1
            message('[',counter,'] ', names(y_logicals)[[counter]] )
         }
         unlist(lapply(as.data.frame(x), function(feature){ main(feature, y_logical) }))
      }))
      m_q_values <- apply(m_p_values, 2, p.adjust, method='bonferroni')

      if(!is.null(max.pvalue)){
         keep_features <- unlist(apply(m_p_values, 2, function(i){
            names(i)[ i<max.pvalue ]
         }), use.names=F)
      } else {
         keep_features <- unlist(apply(m_q_values, 2, function(i){
            names(i)[ i<max.qvalue ]
         }), use.names=F)
      }

      keep_features <- unique(na.exclude(keep_features))
   }

   if(return.new.x){
      return(x[,keep_features,drop=F])
   }
   return(keep_features)
}






