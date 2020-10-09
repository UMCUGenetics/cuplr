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

   if(is.logical(df[,colname.response])){
      y <- df[,colname.response]
      y_ohe <- NA
   } else {
      y <- as.factor(df[,colname.response])
      y_ohe <- oneHotEncode(y, sample.names=rownames(x))
   }

   list(x=x, y=y, y_ohe=y_ohe)
}






