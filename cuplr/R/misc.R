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


