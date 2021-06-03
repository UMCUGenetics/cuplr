#' Execute an extraction function from featureExtractor over a list of files
#'
#' @param file.paths A vector of file paths that will be passed to `func`
#' @param sample.names A vector of sample names
#' @param func A function to execute
#' @param func.args Other args for `func`
#' @param verbose Show progress messages?
#'
#' @return A matrix
#' @export
#'
#' @examples
#' manifest <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/mk_manifest/manifest.txt.gz')
#' if(dir.exists('/Users/lnguyen/')){
#'    counter <- 0
#'    manifest <- as.data.frame(lapply(manifest, function(i){
#'       counter <<- counter + 1
#'       if(counter==1){ return(i) }
#'       paste0('/Users/lnguyen/',i)
#'    }))
#' }
#' 
#' file_paths <- manifest$linx.fusion[1:10]
#' sample_names <- sapply(strsplit(basename(file_paths),'[.]'),`[[`,1)
#' m <- mkFeatureMatrix(
#'    file.paths=file_paths,
#'    sample.names=sample_names,
#'    func=getGeneFusions
#' )
#' 
#' write.table(
#'    m,
#'    gzfile('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/mk_matrices/fusions/m_fusions.txt.gz'),
#'    sep='\t',quote=F
#' )
mkFeatureMatrix <- function(file.paths, sample.names, func, func.args=list(), verbose=T){
   # ## Debugging ------------------------
   # if(F){
   #    manifest <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/mk_manifest/manifest.txt.gz')
   #    if(dir.exists('/Users/lnguyen/')){
   #       counter <- 0
   #       manifest <- as.data.frame(lapply(manifest, function(i){
   #          counter <<- counter + 1
   #          if(counter==1){ return(i) }
   #          paste0('/Users/lnguyen/',i)
   #       }))
   #    }
   #    file.paths=manifest$linx.fusion[1:10]
   #    sample.names=sapply(strsplit(basename(file.paths),'[.]'),`[[`,1)
   #    func=getGeneFusions
   # }
   
   ## Main ------------------------
   if(length(file.paths)!=length(sample.names)){
      stop('`file.paths` and `sample.names` must be the same length')
   }
   
   if(verbose){ message('Applying function to sample:') }
   l <- lapply(1:length(file.paths), function(i){
      #i=1
      if(verbose){ message('[',i,']: ',sample.names[i]) }
      func_args <- c(list(file.paths[i]), func.args)
      do.call(func, func_args)
   })
   
   if(verbose){ message('Merging vectors') }
   m <- do.call(rbind, unname(l))
   rownames(m) <- sample.names
   
   return(m)
}

# mkFeatureMatrix(
#    file.paths=manifest$linx.fusion[1:10],
#    sample.names=sapply(strsplit(basename(file.paths),'[.]'),`[[`,1),
#    func=getGeneFusions
# )

####################################################################################################
#' Fill a matrix with a specified value for missing samples
#'
#' @param m An input matrix
#' @param sample.names A character vector of all samples
#' @param fill.value What value to fill the missing sample rows (default is NA)
#'
#' @return A matrix
#' @export
#'
fillMissingSamplesInMatrix <- function(m, sample.names, fill.value=NA){
   #m=l[[1]]
   #sample.names=all_samples
   
   if(is.null(rownames(m))){ stop('Input matrix must have rownames') }
   if(!all(table(rownames(m) %in% sample.names))){
      warning('Some rownames in the input matrix are not in `sample.names`')
   }
   
   missing_samples <- sample.names[!(sample.names %in% rownames(m))]
   
   if(length(missing_samples)==0){ return(m) }
   
   m_missing <- matrix(
      fill.value, nrow=length(missing_samples), ncol=ncol(m),
      dimnames=list(missing_samples, colnames(m))
   )
   
   return(rbind(m, m_missing))
}

#fillMissingSamplesInMatrix(l[[3]], sample.names=all_samples, fill.value=0)

