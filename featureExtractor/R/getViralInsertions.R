#' Determines the presence of viral inserts in a sample
#'
#' @param viral.inserts Path to the LINX viral inserts txt file, or a dataframe of the file
#'   
#' @return A named logical vector indicating which viruses are present in a sample
#' @export
#'
getViralInsertions <- function(viral.inserts){
   #viral.inserts='/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor//DO48987-from-jar/linx14/DO48987T.linx.viral_inserts.tsv'

   if(is.character(viral.inserts)){
      viral_ins <- read.delim(viral.inserts, check.names=F)
   } else if(is.data.frame(viral.inserts)){
      viral_ins <- viral.inserts
      rm(viral.inserts)
   } else {
      stop('`viral.inserts` must be a path or a dataframe')
   }
   
   virus_group_regex <- c(
      AAV=	'^Adeno-associated virus',
      EBV=	'^Human gammaherpesvirus',
      HBV=	'^(HBV|Hepatitis B)',
      HCV=	'^Hepatitis C virus',
      HIV=	'^Human immunodeficiency virus',
      HPV=	'papillomavirus',
      HSV=	'^Human.*herpesvirus',
      HTLV=	'^Human T-(cell|lymphotropic)',
      MCPyV='^Merkel cell polyomavirus$'
   )
   
   virus_present <- sapply(virus_group_regex, function(i){
      any(grepl(i,viral_ins$VirusName))
   })
   
   return(virus_present)
}


