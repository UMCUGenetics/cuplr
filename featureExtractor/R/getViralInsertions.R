#' Determines the presence of viral inserts in a sample
#'
#' @param linx.viral.inserts.path Path to txt file with the columns: SampleId, 
#' SvId, VirusId, VirusName
#' @param viral.host.ref.path Path to txt file with the columns: refseq_id,
#' virus_name, feature_group, cancer_associated
#' @param return.counts If TRUE, returns the actual number of viral inserts
#' found in a sample
#'   
#' @return An integer vector
#' @export
#'
getViralInsertions <- function(
   linx.viral.inserts.path,
   viral.host.ref.path=VIRAL_HOST_REF,
   return.counts=F, tag.features=T
){
   #linx.viral.inserts.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx/160830_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.viral_inserts.tsv'
   #viral.host.ref.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/svLinxFeatureExtractor/inst/viral_host_ref_ann.txt'

   viral_ref <- read.delim(viral.host.ref.path, check.names=F, stringsAsFactors=F)
   viral_ref <- viral_ref[viral_ref$cancer_associated==1,]

   counts <- sort(unique(viral_ref$feature_group))
   counts <- structure(rep(FALSE,length(counts)), names=counts)

   viral_ins <- read.delim(linx.viral.inserts.path)
   if(nrow(viral_ins)!=0){
      viral_ins$feature_group <- viral_ref[
         match(viral_ins$VirusId, viral_ref$refseq_id), 'feature_group'
      ]
      tab <- table(viral_ins$feature_group)
      counts[ names(tab) ] <- tab
   }
   
   if(return.counts){
      out <- counts
   } else {
      out <- structure(counts>0,names=names(counts))
   }
   
   return(out)
}


