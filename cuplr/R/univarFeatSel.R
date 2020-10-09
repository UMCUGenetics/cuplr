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
#' @param min.effect.size Effect size threshold for keeping features. Default: 3
#' (applies to both +3 and -3 effect sizes)
#' @param min.effect.size.support Minimum number of samples supporting the +ve/-ve effect sizes.
#' Default: 5
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
   max.pvalue=0.01, min.effect.size=3, min.effect.size.support=5, sel.top.n.features=NULL,
   output.type='new.x', verbose=F
){
   if(F){
      max.pvalue=0.01
      min.effect.size=3
      min.effect.size.support=5
      verbose=T

      v.alternative <- rep('greater', ncol(x))
      v.alternative[ grep('(^purple)|(^rmd)',colnames(x)) ] <- 'two.sided'
   }

   if(!is.logical(y)){ stop('`y` must be a logical vector') }

   ## Different test per data type --------------------------------
   pairwiseTest <- function (x, ...) {
      UseMethod("pairwiseTest", x)
   }

   pairwiseTest.numeric <- function(v, alternative='two.sided'){
      #v=df$sigs.snv.SBS15

      a <- v[y]
      b <- v[!y]

      ## Contingency matrix
      amid <- median(a)
      bmid <- median(b)

      apos <- sum(v >  bmid)
      aneg <- sum(v <= bmid)
      bpos <- sum(v >  amid)
      bneg <- sum(v <= amid)

      ## Test
      data.frame(
         apos,aneg,bpos,bneg,
         effect_size=effectSize(apos,aneg,bpos,bneg),
         pvalue=wilcox.test(a, b, alternative=alternative)$p.value,
         alternative
      )
   }

   pairwiseTest.logical <- function(v, alternative='two.sided'){
      #v=df$fusion.TMPRSS2_ERG

      apos <- sum(y & v)
      aneg <- sum(y & !v)

      bpos <- sum(!y & v)
      bneg <- sum(!y & !v)

      data.frame(
         apos,aneg,bpos,bneg,
         effect_size=effectSize(apos,aneg,bpos,bneg),

         ## Fisher test with a modified contingency matrix
         pvalue=fisher.test( matrix(c(apos, aneg+apos, bpos, bneg+bpos),nrow=2), alternative=alternative )$p.value,

         alternative
      )
   }

   pairwiseTest.factor <- function(v, alternative='two.sided'){
      #v=x$purple.gender

      neg_category <- levels(v)[1]

      apos <- sum(y & v!=neg_category)
      aneg <- sum(y & v==neg_category)

      bpos <- sum(!y & v!=neg_category)
      bneg <- sum(!y & v==neg_category)

      ## Fisher test with a modified contingency matrix
      data.frame(
         apos,aneg,bpos,bneg,
         effect_size=effectSize(apos,aneg,bpos,bneg),

         ## Fisher test with a modified contingency matrix
         pvalue=fisher.test( matrix(c(apos, aneg+apos, bpos, bneg+bpos),nrow=2), alternative=alternative )$p.value,

         alternative
      )
   }

   #pairwiseTest(df$sigs.snv.SBS15)
   #pairwiseTest(df$fusion.TMPRSS2_ERG)
   #pairwiseTest(x$purple.gender)

   ## Main --------------------------------

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

   counter <- 0
   tests <- do.call(rbind, lapply(x, function(i){
      counter <<- counter + 1
      if(verbose){ message('[',counter,'] ', colnames(x)[[counter]] ) }
      pairwiseTest(i, alternative=v.alternative[counter])
   }))
   tests <- tests[order(tests$pvalue),]

   ## Select features
   tests_ss <- subset(
      tests,
      pvalue < max.pvalue
   )

   if(!is.null(min.effect.size) & !is.null(min.effect.size.support)){
      tests_ss <- subset(
         tests_ss,
            (alternative %in% c('two.sided','greater') & effect_size >=  min.effect.size & apos >= min.effect.size.support) |
            (alternative %in% c('two.sided','less')    & effect_size <= -min.effect.size & bpos >= min.effect.size.support)
      )
   }

   keep_features <- rownames(tests_ss)

   if(!is.null(sel.top.n.features)){
      keep_features <- keep_features[1:sel.top.n.features]
      keep_features <- na.exclude(keep_features)
   }

   if(output.type=='raw'){
      return(tests)
   } else if(output.type=='features'){
      return(keep_features)
   }
   return(x[,keep_features,drop=F])
}


# univarFeatSel <- function(
#    x, y,
#    max.qvalue=0.01, max.pvalue=NULL, sel.top.n.features=NULL,
#    return.new.x=T, verbose=F
# ){
#    # colname.response='response'
#    # x <- df[,colnames(df)!=colname.response]
#    # y <- as.factor(df[,colname.response])
#    # y <- y=='Prostate'
#
#    if( !(is.logical(y) | is.factor(y)) ){
#       stop('`y` must be a logical or factor')
#    }
#    #if(is.logical(y)){ y <- factor(y,c('TRUE','FALSE')) }
#
#    main <- function(v, y.logical){
#       #y.logical=y
#
#       ## Numeric data: wilcox text
#       if(is.numeric(v)){
#          #v=x$viral_ins.Hepatitis_C_virus
#          #v=x$rmd.14q_107
#          v_split <- split(v, y.logical)
#          wilcox.test(v_split[['TRUE']], v_split[['FALSE']])$p.value
#
#          ## Categorical data: fisher test
#       } else {
#          #v=ifelse(x$purple.gender,'male','female')
#          #v=unname(m[,'AR'])
#          #v=x$gene_def.VHL
#          #y=metadata[ match(rownames(m), metadata$sample),'cancer_type' ]
#          #y.logical <- y=='Prostate'
#
#          neg_category <- levels(as.factor(v))[1]
#          #v=rep('0;none',length(v))
#
#          fisher.test(
#             matrix(
#                c(
#                   sum(v!=neg_category & y.logical), sum(y.logical),
#                   sum(v!=neg_category & !y.logical), length(y.logical)
#                ),
#                nrow=2
#             )
#          )$p.value
#       }
#    }
#
#    ## Binary classification
#    if(is.logical(y)){
#       if(verbose){ counter <- 0 }
#       p_values <- unlist(lapply(as.data.frame(x), function(i){
#          if(verbose){
#             counter <<- counter + 1
#             message('[',counter,'] ', colnames(x)[[counter]] )
#          }
#          main(i, y)
#       }))
#       p_values <- sort(p_values)
#       q_values <- p.adjust(p_values, method='bonferroni')
#
#       if(!is.null(max.pvalue)){
#          keep_features <- names(p_values)[ p_values < max.pvalue ]
#       } else {
#          keep_features <- names(q_values)[ q_values < max.qvalue ]
#       }
#
#       if(!is.null(sel.top.n.features)){ keep_features <- keep_features[1:sel.top.n.features] }
#       keep_features <- na.exclude(keep_features)
#
#       ## Multiclass classification
#    } else {
#       y_logicals <- lapply(levels(y), function(i){ y==i })
#       names(y_logicals) <- levels(y)
#
#       if(verbose){ counter <- 0 }
#       m_p_values <- do.call(cbind, lapply(y_logicals, function(y_logical){
#          if(verbose){
#             counter <<- counter + 1
#             message('[',counter,'] ', names(y_logicals)[[counter]] )
#          }
#          unlist(lapply(as.data.frame(x), function(feature){ main(feature, y_logical) }))
#       }))
#       m_q_values <- apply(m_p_values, 2, p.adjust, method='bonferroni')
#
#       if(!is.null(max.pvalue)){
#          keep_features <- unlist(apply(m_p_values, 2, function(i){
#             names(i)[ i<max.pvalue ]
#          }), use.names=F)
#       } else {
#          keep_features <- unlist(apply(m_q_values, 2, function(i){
#             names(i)[ i<max.qvalue ]
#          }), use.names=F)
#       }
#
#       keep_features <- unique(na.exclude(keep_features))
#    }
#
#    if(return.new.x){
#       return(x[,keep_features,drop=F])
#    }
#    return(keep_features)
# }





