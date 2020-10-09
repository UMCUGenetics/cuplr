#' Count number of insertions of repeat elements
#'
#' @param vcf.file Path to vcf file
#' @param rep.elem.whitelist.path Path to txt file with list of selected repeat elements
#'
#' @return An integer vector
#' @export
#'
countRepElemInsertions <- function(
   vcf.file, rep.elem.whitelist.path=REP_ELEM_WHITELIST, 
   tag.features=T
){
   #vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/160704_HMFregCPCT_FR12244543_FR12244602_CPCT02030245/CPCT02030245T.purple.sv.ann.vcf.gz'
   #vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/171223_HMFregCPCT_FR15580187_FR15579534_CPCT02010704/CPCT02010704T.purple.sv.ann.vcf.gz'
   
   vcf <- mutSigExtractor::variantsFromVcf(
      vcf.file, vcf.filter='PASS', 
      vcf.fields=c('CHROM','POS','REF','ALT','FILTER','INFO'),
      keep.chroms=c(1:22,'X')
   )
   
   rep_elem_whitelist <- read.table(rep.elem.whitelist.path, header=F, stringsAsFactors=F)[,1]
   counts <- structure(rep(0,length(rep_elem_whitelist)), names=rep_elem_whitelist)
   
   if(!is.data.frame(vcf)){ return(counts) }
   
   rep_elements <- mutSigExtractor::getInfoValues(vcf$info, 'INSRMRC')[,1]
   
   tab <- table(rep_elements)
   tab <- tab[names(tab) %in% rep_elem_whitelist]
   
   counts[names(tab)] <- tab
   
   return(counts)
}
