#' Extract SV features from PURPLE and LINX output
#'
#' @param purple.cnv.path Path to purple cnv file
#' @param purple.purity.path Path to purple.purity.tsv file
#' @param linx.fusion.path Path to LINX txt file containing fusion data
#' @param linx.viral.inserts.path Path to txt file with the columns: SampleId, 
#' SvId, VirusId, VirusName
#' @param driver.catalog.path Path to LINX driver catalog txt file
#' @param sample.name See out.path
#' @param out.path If provided, will write output to path as a 1 column dataframe with the 
#' sample.name as the column name and feature names as the row names
#' @param verbose Show messages?
#'
#' @return A 1-row data.frame
#' @export
#'
extractSvFeatures <- function(
   purple.cnv.path, purple.purity.path,
   linx.fusion.path, linx.viral.inserts.path, driver.catalog.path,
   out.path=NULL, sample.name=NULL, verbose=1
){
   ## Debugging --------------------------------
   if(F){
      manifest <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/mk_manifest/manifest.txt.gz')
      sample_name <- 'CPCT02020731T'
      
      counter <- 0
      manifest <- as.data.frame(lapply(manifest, function(i){
         counter <<- counter + 1
         if(counter==1){ return(i) }
         paste0('/Users/lnguyen/',i)
      }))
      
      purple.cnv.path=manifest[manifest$sample==sample_name,'purple.cnv']
      purple.purity.path=manifest[manifest$sample==sample_name,'purple.purity']
      linx.fusion.path=manifest[manifest$sample==sample_name,'linx.fusion']
      linx.viral.inserts.path=manifest[manifest$sample==sample_name,'linx.viral_inserts']
      driver.catalog.path=manifest[manifest$sample==sample_name,'driver.catalog']
      
      verbose=2
      
      devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/')
   }
   
   ##--------------------------------
   features <- list()
   
   if(verbose){ message('> Getting PURPLE gender and whole genome duplication status...') }
   features$purple_purity <- getPurplePurityData(purple.purity.path)
   
   ##--------------------------------
   if(verbose){ message('> Calculating chrom arm ploidies...') }
   features$ploidy <- calcChromArmPloidies(purple.cnv.path, verbose=(verbose==2))
   
   ##--------------------------------
   if(verbose){ message('> Identifying gene fusions...') }
   features$fusions <- getGeneFusions(linx.fusion.path)
   
   ##--------------------------------
   if(verbose){ message('> Identifying viral inserts...') }
   features$viral_inserts <- getViralInsertions(linx.viral.inserts.path)
   
   ##--------------------------------
   if(verbose){ message('> Getting gene amplification levels...') }
   features$gene_amp <- getGeneAmpLevels(driver.catalog.path)
   
   ##--------------------------------
   if(verbose){ message('> Calculating gene mut scores...') }
   features$gene_mut <- calcGeneMutScores(driver.catalog.path)
   
   ##--------------------------------
   if(verbose){ message('> Making final feature vector...') }
   
   out <- do.call(c, unname(features))

   if(is.null(out.path)){
      return(out)
   } else {
      sample_name <- if(!is.null(sample.name)){ sample.name } else { basename(purple.cnv.path) }
      write.table(
         matrix( out, dimnames=list(names(out), sample_name) ), 
         out.path, sep='\t', row.names=F, quote=F
      )
   }
   
   # out <- lapply(features, function(i){
   #    if(is.vector(i)){ as.data.frame(rbind(i)) }
   #    else { return(i) }
   # })
   # 
   # out <- do.call(cbind, unname(out))
   # rownames(out) <- NULL
   # 
   # if(is.null(out.path)){
   #    return(out)
   # } else {
   #    write.table(out, out.path, sep='\t', row.names=F, quote=F)
   # }
}


####################################################################################################
#' Execute an extraction function from hmfSvFeatureExtractor over a list of files
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








