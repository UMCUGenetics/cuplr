#' Get gender and whole genome duplication status
#'
#' @param purple.purity.path Path to purple.purity.tsv file
#'
#' @return A dataframe
#' @export
#'
getPurplePurityData <- function(purple.purity.path){
   #purple.purity.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/160704_HMFregXXXXXXXX/XXXXXXXX.purple.purity.tsv'
   df <- read.delim(purple.purity.path, check.names=F, stringsAsFactors=F)

   c(
      gender = df$gender=='FEMALE', ## MALE_KLINEFELTER becomes male
      genome_ploidy = df$ploidy,
      diploid_proportion = df$diploidProportion,
      has_wgd = as.logical(df$wholeGenomeDuplication)
   )
}

