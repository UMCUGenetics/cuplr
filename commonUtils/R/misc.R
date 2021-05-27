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

####################################################################################################
#' Return first rows and columns of a matrix or dataframe
#'
#' @param m A matrix or dataframe
#' @param rows Number of rows
#' @param cols Number of cols
#'
#' @return The subsetted matrix or dataframe
#' @export
#'
head2 <- function(m, rows=10, cols=10){

   if(nrow(m)<rows){ rows <- nrow(m) }
   if(ncol(m)<cols){ cols <- ncol(m) }

   m[1:rows,1:cols]
}

####################################################################################################
#' Convert all columns of a dataframe into factors
#'
#' @param df A dataframe
#' @param levels A character vector of factor levels
#' @param ... Arguments that can be passed to factor()
#'
#' @return A factor dataframe
#' @export
#'
as.factor.data.frame <- function(df, levels, ...){
   #df=features$gene_def
   #levels=c("none","mut","mut,mut","loh,mut","loh_arm,mut","loh_chrom,mut","deep_deletion")

   out <- as.data.frame(lapply(df, function(i){set.seed
      factor(i, levels, ...)
   }))
   dimnames(out) <- dimnames(df)
   return(out)
}


####################################################################################################
#' Cache remote files and read as a dataframe
#'
#' @param remote.path Remote path to txt or rds file
#' @param local.path Local path to copy to
#' @param cache.dir Path to the cache directory. If `local.path` is not specified, files will be
#' cached here.
#' @param overwrite Redownload data from remote?
#'
#' @return A dataframe
#' @export
#'
cacheAndReadData <- function(
   remote.path, local.path=NULL,
   cache.dir=path.expand('~/Documents/R_cache/'),
   overwrite=F
){

   #remote.path=paste0(base_dir,'/misc/processed/Chromatin_modifiers/scripts/annotate_svs/vis_sv_data_compact.txt.gz')

   ## Init --------------------------------
   ext <- c('.rds','.txt','.txt.gz','.csv','.csv.gz')
   regex <- gsub('[.]','[.]',ext)
   regex <- paste0(regex,'$')

   is_valid_ext <- sapply(regex, function(i){ grepl(i, remote.path) })
   if(all(!is_valid_ext)){
      stop('File extension must be one of the following:\n  ', paste(ext, collapse=', '))
   }

   ext <- ext[ which(is_valid_ext) ]
   regex <- regex[ which(is_valid_ext) ]

   set.seed(nchar(remote.path))

   ## Copy --------------------------------
   if(dir.exists('/hpc/')){
      local.path <- remote.path
   } else {
      if(is.null(local.path)){
         local.path <- paste0(
            cache.dir,
            sub(regex,'',basename(remote.path)),
            '.',paste(sample(letters, 8), collapse=''),ext
         )
      }

      local.path <- gsub('/+', '/', local.path)
      if(!file.exists(local.path) | overwrite){
         if(!file.exists(local.path)){
            message('Making local copy: ', local.path)
         } else {
            message('Updating local copy: ', local.path)
         }

         system(sprintf(
            'rsync -a %s %s',
            remote.path,
            local.path
         ))
      }
   }

   ## Read --------------------------------
   message('Reading local copy: ', local.path)
   if(ext=='.rds'){
      readRDS(local.path)
   } else if(grepl('txt',ext)){
      read.delim(local.path, check.names=F)
   } else {
      read.csv(local.path, check.names=F)
   }
}


