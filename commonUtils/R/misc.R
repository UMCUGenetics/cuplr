#' Assign clusters to a vector of values based on run length encoding
#'
#' @param x An integer vector
#'
#' @return
#' @export
#'
rle2Clusters <- function(x){
   rle_out <- rle(x)
   unlist(lapply(1:length(rle_out$lengths), function(i){
      rep(i, rle_out$lengths[i])
   }))
}

####################################################################################################
#' Combine feature columns by summing them together
#'
#' @param x A matrix or dataframe
#' @param target.features A character vector of column (feature) names
#' @param regex Instead of target.features, a regex can be specified
#' @param target.name A name to assign the new feature. If unspecified, will use the name of the
#' first old feature
#'
#' @return The original matrix or dataframe with the indicated features combined
#' @export
#'
combineFeatures <- function(x, target.features=NULL, regex=NULL, target.name=NULL){

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

####################################################################################################
#' Insert a row into a dataframe
#'
#' @param df A dataframe
#' @param new.row The row as a vector or dataframe to insert
#' @param row.num Row index to insert the row. Rows below this are shifted down.
#' @param offset Insert the row at +offset rows after row.num
#'
#' @return A dataframe
#' @export
#'
insertRow <- function(df, row.num, new.row=NULL, offset=1) {
   df[seq(row.num+1,nrow(df)+1),] <- df[seq(row.num,nrow(df)),]

   if(is.null(new.row)){
      new.row <- df[row.num,]
   }

   df[row.num+offset,] <- new.row

   return(df)
}

####################################################################################################
#' Write tsv file
#'
#' @param x The object to be written, preferably a matrix or data frame. If not, it is attempted to
#' coerce x to a data frame.
#' @param file Path to the output file. If file ends with .gz, a gzip compressed file will be written
#' @param ... Arguments that can be passed to write.table()
#'
#' @return
#' @export
#'
write.tsv <- function(x, file, ...){
   write.table(
      x,
      if(grepl('[.]gz$',file)){ gzfile(file) } else { file },
      sep='\t', quote=F, row.names=F,
      ...
   )
}

####################################################################################################
#' Insert column(s) after a column of an existing dataframe
#'
#' @param df Input dataframe
#' @param v A vector or dataframe
#' @param after Column index or name
#' @param colname Set inserted column name. Has no effect if v is not a vector
#'
#' @return A dataframe with the inserted column(s)
#' @export
#'
insColAfter <- function(df, v, after, colname=NULL){
   if(is.character(after)){
      after <- which(colnames(df)==after)
   }

   df_l <- df[1:after]
   df_r <- df[(after+1):ncol(df)]

   df_new <- cbind(df_l, v, df_r)

   if(!is.null(colname)){
      colnames(df_new)[(after+1)] <- colname
   }

   return(df_new)
}








