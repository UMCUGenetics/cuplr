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


   ## Force factor levels from training set onto new data
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












