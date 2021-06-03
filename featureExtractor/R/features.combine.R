#' Combine feature columns by summing them together
#'
#' @rdname combineFeatures
#'
#' @param x An integer/numeric vector, matrix or dataframe
#' @param target.features A character vector of column (feature) names
#' @param regex Instead of target.features, a regex can be specified
#' @param target.name A name to assign the new feature. If unspecified, will use the name of the
#' first old feature
#'
#' @return The original data with the indicated features combined
#' @export
#'
combineFeatures <- function (x, ...) {
   UseMethod("combineFeatures", x)
}

#' @rdname combineFeatures
#' @method combineFeatures matrix
combineFeatures.matrix <- function(x, target.features=NULL, regex=NULL, target.name=NULL){

   if(!is.null(regex)){
      target.features <- grep(regex,colnames(x),value=T)
      if(length(target.features)==0){ stop('No features match the regex pattern') }
   }

   target_col <- target.features[1]

   x[,target_col] <- rowSums(x[,target.features])

   rm_cols <- target.features[target.features!=target_col]
   x <- x[,!(colnames(x) %in% rm_cols)]

   if(!is.null(target.name)){
      colnames(x)[colnames(x)==target_col] <- target.name
   }

   return(x)
}

#' @rdname combineFeatures
#' @method combineFeatures dataframe
combineFeatures.data.frame <- combineFeatures.matrix

#' @rdname combineFeatures
#' @method combineFeatures numeric
combineFeatures.numeric <- function(x, target.features=NULL, regex=NULL, target.name=NULL){

   if(!is.null(regex)){
      target.features <- grep(regex,names(x),value=T)
      if(length(target.features)==0){ stop('No features match the regex pattern') }
   }

   target_index <- target.features[1]

   x[target_index] <- sum(x[target.features])

   rm_cols <- target.features[target.features!=target_index]
   x <- x[!(names(x) %in% rm_cols)]

   if(!is.null(target.name)){
      names(x)[names(x)==target_index] <- target.name
   }

   return(x)
}

#' @rdname combineFeatures
#' @method combineFeatures integer
combineFeatures.integer <- combineFeatures.numeric
