#' Determine the presence of gene fusions
#'
#' @param linx.fusion.path Path to LINX txt file containing fusion data
#' @param gene.fusion.whitelist.path Path to txt file with list of gene fusions
#' @param tag.features If TRUE, will tag the value names with the feature type
#'
#' @return A (0,1) vector
#' @export
#' 
getGeneFusions <- function(
   linx.fusion.path, gene.fusion.whitelist.path=GENE_FUSION_WHITELIST,
   tag.features=T
){
   #linx.fusion.path='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180502_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.fusion.tsv'
   linx_fusions <- read.delim(linx.fusion.path, check.names=F, stringsAsFactors=F)
   gene_fusion_whitelist <- read.table(gene.fusion.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   out <- structure(
      as.integer(gene_fusion_whitelist %in% linx_fusions$Name),
      names=gene_fusion_whitelist
   )
   
   if(tag.features){
      names(out) <- paste0('fusion.',names(out))
   }
   
   return(out)
}

