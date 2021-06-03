#' Group a vector of feature names by their feature tags
#'
#' @param x A character vector in the form 'feature_tag.feature_name'
#' @param rm.tags Remove feature group tags?
#'
#' @return A list of feature vectors in for each tag
#' @export
#'
groupFeaturesByTag <- function(x, rm.tags=F){
   feature_tags <- gsub('[.].*$','',x)

   out <- split(x, factor(feature_tags, unique(feature_tags)))
   if(rm.tags){
      out <- lapply(out, function(i){ gsub('^\\w+[.]','',i) })
   }

   return(out)
}

####################################################################################################
#' Split a dataframe into feature groups
#'
#' @rdname splitFeaturesByGroup
#'
#' @param df A dataframe where colnames are in the form
#' 'feature_tag.feature_name'
#' @param rm.tags Remove feature group tags?
#'
#' @return A list of dataframes
#' @export
#'
splitFeaturesByGroup <- function (x, ...) {
   UseMethod("splitFeaturesByGroup", x)
}

#' @rdname splitFeaturesByGroup
#' @method splitFeaturesByGroup default
#' @export
splitFeaturesByGroup.default <- function(x, rm.tags=F){
   feature_groups <- groupFeaturesByTag(names(x))
   l <- lapply(feature_groups, function(i){ x[i] })
   if(rm.tags){
      l <- lapply(l, function(i){
         names(i) <- gsub('^\\w+[.]','',names(i))
         return(i)
      })
   }
   return(l)
}

#' @rdname splitFeaturesByGroup
#' @method splitFeaturesByGroup data.frame
#' @export
splitFeaturesByGroup.data.frame <- function(x, rm.tags=F){
   feature_groups <- groupFeaturesByTag(colnames(x))
   l <- lapply(feature_groups, function(i){ x[,i,drop=F] })
   if(rm.tags){
      l <- lapply(l, function(i){
         colnames(i) <- gsub('^\\w+[.]','',colnames(i))
         return(i)
      })
   }
   return(l)
}


