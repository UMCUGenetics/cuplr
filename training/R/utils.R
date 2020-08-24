################################################################################
#' Preserve order of dataframe values
#'
#' @description Converts character columns into factors to preserve the order of their values
#' @param df A dataframe
#'
#' @return A data frame with characters converted into factors
#' @export
#'
forceDfOrder <- function(df){
   if(!is.data.frame(df)){ stop('Input must be a dataframe') }
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

################################################################################
#' Prepares a feature dataframe for multiclass classifier training
#'
#' @description Separates a dataframe containing features and response variable. The response
#' variable is also one hot encoded
#'
#' @param df A dataframe
#' @param colname.response Name of the response column
#'
#' @return A list containing the feature matrix/dataframe, response vector, and response
#' one-hot encode matrix
#' @export
#'
dfToFeaturesAndResponse <- function(df, colname.response='response'){
   x <- df[,colnames(df)!=colname.response]
   #x <- as.matrix(x)

   y <- as.factor(df[,colname.response])
   y_ohe <- oneHotEncode(y, sample.names=rownames(x))

   list(x=x, y=y, y_ohe=y_ohe)
}

################################################################################
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
      if(is.numeric(v)){
         #v=x$viral_ins.Hepatitis_C_virus
         #v=x$rmd.14q_107
         v_split <- split(v, y.logical)
         wilcox.test(v_split[['TRUE']], v_split[['FALSE']])$p.value
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
         keep_features <- names(p_values)[ p_values < max.qvalue ]
      } else {
         keep_features <- names(q_values)[ q_values < max.qvalue ]
      }

      if(!is.null(sel.top.n.features)){ keep_features <- keep_features[1:sel.top.n.features] }
      keep_features <- na.exclude(keep_features)

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

      keep_features <- unlist(apply(m_q_values, 2, function(i){
         names(i)[ i<max.qvalue ]
      }), use.names=F)

      keep_features <- unique(na.exclude(keep_features))
   }

   if(return.new.x){
      return(x[,keep_features,drop=F])
   }
   return(keep_features)
}



################################################################################
#' Predict method for random forest ensemble
#'
#' @param object An object of class randomForestEnsemble
#' @param newdata A data frame or matrix containing new data. (Note: If not given, the out-of-bag
#' prediction in object is returned)
#' @param type 'response', 'prob' or 'votes', indicating the type of output: predicted values,
#' matrix of class probabilities, or matrix of vote counts. 'class' is allowed, but automatically
#' converted to "response", for backward compatibility.
#'
#' @return 'response': predicted classes (the classes with majority vote). 'prob' matrix of class
#' probabilities (one column for each class and one row for each input). 'vote'
#' matrix of vote counts (one column for each class and one row for each new input); either in raw
#' counts or in fractions (if norm.votes=TRUE).
#' @export
#'
predict.randomForestEnsemble <- function(object, newdata, type='response', verbose=F){
   # object=model
   # newdata=test_data$x
   # type='response'

   categorical_feature_names <- names(object[[1]]$categorical_lvls)
   x <- as.data.frame(lapply(colnames(newdata), function(i){
      #i='gene_def.AR'
      if(!(i %in% categorical_feature_names)){
         return(newdata[,i])
      } else{
         factor(newdata[,i], levels=object[[1]]$categorical_lvls[[i]])
      }
   }))
   rownames(x) <- rownames(newdata)
   colnames(x) <- colnames(newdata)

   if(verbose){ counter <- 0 }
   m <- do.call(cbind, lapply(object, function(i){
      #i=object$Breast
      if(verbose){
         counter <<- counter + 1
         message('[RF ',counter,'/',length(object),']: ', names(object)[counter])
      }
      randomForest:::predict.randomForest(i, x, type='prob')[,1]
   }))

   if(type=='prob'){
      m
   } else {
      factor( colnames(m)[ max.col(m) ], levels=colnames(m) )
   }
}












