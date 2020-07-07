#' Group a vector of feature names by their feature tags
#' 
#' @param x A character vector in the form 'feature_tag.feature_name'
#' @param rm.tags Remove feature group tags?
#'
#' @return A list of feature vectors in for each tag
#' @export
#'
groupFeaturesByTag <- function(x, rm.tags=F){
   #df <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/extract_sv_features/svlinx_features.txt.gz')
   #x=colnames(df)
   feature_tags <- regmatches(x, regexpr('^\\w+[.]', x))
   feature_tags <- gsub('[.]$','',feature_tags)
   
   out <- split(x, factor(feature_tags, unique(feature_tags)))
   if(rm.tags){
      out <- lapply(out, function(i){ gsub('^\\w+[.]','',i) })
   }
   
   return(out)
}

#' Split a dataframe into feature groups
#' 
#' @param df A dataframe where colnames are in the form 
#' 'feature_tag.feature_name'
#' @param rm.tags Remove feature group tags?
#'
#' @return A list of dataframes
#' @export
#'
splitDfByFeatureGroup <- function(df, rm.tags=F){
   feature_groups <- groupFeaturesByTag(colnames(df))
   l <- lapply(feature_groups, function(i){ df[,i] })
   if(rm.tags){
      l <- lapply(l, function(i){
         colnames(i) <- gsub('^\\w+[.]','',colnames(i))
         return(i)
      })
   }
   return(l)
}
