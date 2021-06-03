#' Determine the presence of gene mono/biallelic hits, deep deletions, and amplifications
#'
#' @param linx.fusions Path to the LINX driver catalog txt file, or a dataframe of the file
#' @param whitelist.path Path to the gene whitelist txt file
#' @param min.dnds.likelihood Min DNDS likelihood for to consider a monoallelic hit as impactful
#' @param na.output If TRUE, and if `linx.fusions` is NA, then a vector of FALSE values will be 
#' returned
#' 
#' @return A named logical vector indicating which gene driver events are present in a sample
#' @export
#' 
getGeneDriverEvents <- function(
   linx.drivers, whitelist.path=GENE_DRIVER_WHITELIST, 
   min.dnds.likelihood=0.9, na.output=T
){
   
   # if(F){
   #    linx.drivers='/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/DO218121-from-jar/linx14/DO218121T.linx.driver.catalog.tsv'
   #    #linx.drivers='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/gene_drivers/05_HMF_PCAWG_full/driver.catalog.txt.gz'
   #    whitelist.path=GENE_DRIVER_WHITELIST
   #    min.dnds.likelihood=0.9
   # }
   
   ## Init --------------------------------
   ## Initialize output
   driver_whitelist <- read.table(whitelist.path, header=F, stringsAsFactors=F)[,1]
   event_types <- c('amp','deep_del','biall', 'monoall')
   
   driver_event_exists <- unlist(lapply(driver_whitelist, function(i){ paste0(i,'.',event_types) }))
   driver_event_exists <- structure(rep(FALSE, length(driver_event_exists)), names=driver_event_exists)

   ## Load input
   if(is.na(linx.drivers) & na.output){ 
      warning('`linx.drivers` is NA. Returning a vector of FALSE values')
      return(driver_event_exists) 
   }
   
   if(is.character(linx.drivers)){
      drivers <- read.delim(linx.drivers, check.names=F)
   } else if(is.data.frame(linx.drivers)){
      drivers <- linx.drivers
      #rm(linx.drivers)
   } else {
      stop('`linx.drivers` must be a path or a dataframe')
   }
   
   drivers <- subset(drivers, gene %in% driver_whitelist)
   
   if(nrow(drivers)==0){ return(driver_event_exists) }
   
   ## Main --------------------------------
   #df <- drivers
   #df <- df[,c('sample', 'gene', 'driver', 'category', 'likelihoodMethod', 'driverLikelihood', 'biallelic')]
   drivers$biallelic <- as.logical(drivers$biallelic)
   
   drivers$event_type <- 'none'
   drivers$event_type <- with(drivers,{
      
      event_type[driver=='MUTATION' & !biallelic] <- 'monoall'
      event_type[biallelic & likelihoodMethod!='DEL'] <- 'biall'
      event_type[biallelic & likelihoodMethod=='DEL'] <- 'deep_del'
      event_type[likelihoodMethod=='AMP'] <- 'amp'
      
      factor(event_type, c('none','monoall','biall','deep_del','amp'))
   })
   
   driver_events <- paste0(drivers$gene,'.',drivers$event_type)
   driver_events <- driver_events[driver_events %in% names(driver_event_exists)]
   
   driver_event_exists[driver_events] <- TRUE
   
   return(driver_event_exists)
}