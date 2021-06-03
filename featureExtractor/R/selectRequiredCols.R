#' Check column names of an input table
#'
#' @param df An input dataframe
#' @param required.cols A character vector indicating the require col names
#' @param sel.cols A named character vector with the names. Each name corresponds to the names of
#' the `required.cols`. Each value corresponds to the the column names in the txt file. If
#' `sel.cols` is NULL, this will automatically be generated from required.cols
#'
#' @description This function checks whether all `names(sel.cols)` are in `required.cols`. If this
#' passes, the function will check whether all `sel.cols` values are in `colnames(df)`. If this,
#' passes, a dataframe will be returned with the selected columns
#'
#' @return A subset dataframe if all checks pass
#' @export
#'
selectRequiredCols <- function(df, required.cols, sel.cols=NULL){
   #required.cols=c('ResolvedType','ClusterId','PosStart','PosEnd')

   if(!is.data.frame(df) & !is.matrix(df)){
      stop("`df` must be a dataframe or matrix")
   }

   if(is.null(sel.cols) || is.na(sel.cols)){
      sel.cols <- structure(required.cols, names=required.cols)
   }
   if(length(names(sel.cols))==0){
      stop("`sel.cols` must be a named character vector")
   }

   missing_names <- required.cols[ !(required.cols %in% names(sel.cols)) ]
   if(length(missing_names)>0){
      stop("`sel.cols` is missing the names: ", paste(missing_names, collapse=', '))
   }

   missing_cols <- sel.cols[!(sel.cols %in% colnames(df))]
   if(length(missing_cols)>0){
      stop("Some columns are missing in the input table: ",paste(missing_cols, collapse=', '))
   }

   df <- df[,sel.cols,drop=F]
   colnames(df) <- names(sel.cols)
   return(df)
}

