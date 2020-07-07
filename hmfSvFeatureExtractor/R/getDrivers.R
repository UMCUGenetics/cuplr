#' Get minCopyNumber for set of genes
#'
#' @param driver.catalog.path Path to LINX driver catalog txt file
#' @param df A dataframe of the driver catalog txt file
#' @param gene.whitelist.path Path to txt file with list of selected genes
#' @param tag.features If TRUE, will tag the value names with the feature type
#'
#' @return A numeric vector
#' @export
#'
getGeneAmpLevels <- function(
   driver.catalog.path=NULL, df=NULL, gene.whitelist.path=GENE_AMP_WHITELIST, tag.features=T
){
   #driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv'
   #gene.whitelist.path=GENE_AMP_WHITELIST
   
   gene_whitelist <- read.table(gene.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   if(is.null(df)){
      df <- read.delim(driver.catalog.path)
   }
   
   df <- df[df$driver=='AMP' & df$gene %in% gene_whitelist,]
   df <- df[!duplicated(df$gene),c('gene','minCopyNumber'),]
   
   out <- structure(rep(0, length(gene_whitelist)), names=gene_whitelist)
   if(nrow(df)!=0){
      out[df$gene] <- round(df$minCopyNumber, 2)
   }
   
   if(tag.features){
      names(out) <- paste0('gene_amp.',names(out))
   }
   
   return(out)
}


####################################################################################################
#' Calculate score for LINX driver mutation types
#'
#' @param df A dataframe of the driver catalog txt file
#'
#' @return An integer vector
#' @export
#'
#' @examples
#' df <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/select_drivers//driver.catalog.merged.txt.gz')
#' scoreDriverCatalog(df)
scoreDriverCatalog <- function(df){
   #driver.catalog.path=paste0(base_dir,'/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/select_drivers//driver.catalog.merged.txt.gz')
   #df <- read.delim(driver.catalog.path)

   with(df,{
      out <- rep(0, nrow(df))
      out[ driver=='AMP'] <- -1
      
      ## Monoallelic loss (hotspot/DNDS)
      out[ biallelic=='false' & driver=='MUTATION' & likelihoodMethod != 'BIALLELIC'] <- 1
      
      ## Biallelic loss (LOH + hotspot/DNDS)
      out[ biallelic=='true' & driver=='MUTATION' & likelihoodMethod != 'BIALLELIC' ] <- 2
      
      ## Biallelic loss (LOH + known pathogenic variant)
      out[ biallelic=='true' & driver=='MUTATION' & likelihoodMethod == 'BIALLELIC' ] <- 3
      
      ## Biallelic loss (2x SV ?)
      out[ driver=='HOM_DISRUPTION' ] <- 4
      
      ## Deep deletion
      out[ driver=='DEL' ] <- 5
      
      return(out)
   })
}

#' Score mutations for each gene in whitelist
#'
#' @param driver.catalog.path Path to LINX driver catalog txt file
#' @param df A dataframe of the driver catalog txt file
#' @param gene.whitelist.path Path to txt file with list of selected genes
#' @param tag.features If TRUE, will tag the value names with the feature type
#'
#' @return An integer vector
#' @export
#'
calcGeneMutScores <- function(
   driver.catalog.path=NULL, df=NULL, gene.whitelist.path=GENE_MUT_WHITELIST, tag.features=T
){
   #driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv'
   #gene.whitelist.path=GENE_MUT_WHITELIST
   
   gene_whitelist <- read.table(gene.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   if(is.null(df)){
      df <- read.delim(driver.catalog.path)
   }
   
   df <- df[df$driver!='AMP' & df$gene %in% gene_whitelist,]
   df <- df[!duplicated(df$gene),]
   df$score <- scoreDriverCatalog(df)
   
   out <- structure(rep(0, length(gene_whitelist)), names=gene_whitelist)
   if(nrow(df)!=0){
      out[df$gene] <- df$score
   }
   
   if(tag.features){
      names(out) <- paste0('gene_mut.',names(out))
   }
   
   return(out)
}

####################################################################################################
#' Get mutation status from a whitelist of genes
#'
#' @param driver.catalog.path Path to LINX driver catalog txt file
#' @param df A dataframe of the driver catalog txt file
#' @param gene.whitelist.path Path to txt file with list of selected genes
#' @param mode Can be 'biallelic', 'monoallel', or 'amp'
#' @param tag.features If TRUE, will tag the value names with the feature type
#'
#' @return A (1,0 vector)
#' @export
#'
getGeneMutStatus <- function(
   driver.catalog.path=NULL, df=NULL, gene.whitelist.path, mode,
   tag.features=T
){
   #driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv'
   #gene.whitelist.path=GENE_AMP_WHITELIST
   #driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv'
   #gene.whitelist.path=GENE_BIALLEL_WHITELIST
   
   gene_whitelist <- read.table(gene.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   if(is.null(df)){ df <- read.delim(driver.catalog.path) }
   
   df <- df[df$gene %in% gene_whitelist,]
   df$score <- scoreDriverCatalog(df)
   df <- df[order(df$score, decreasing=T),]
   df <- df[!duplicated(df$gene),]
   
   if(mode=='biall'){
      df <- df[df$score>=2,]
   } else if(mode=='monoall'){
      df <- df[df$score==1,]
   } else if(mode=='amp'){
      df <- df[df$score == -1,]
   } else {
      stop("`mode` must be 'biall', 'monoall', or 'amp'")
   }
   
   out <- structure(rep(0, length(gene_whitelist)), names=gene_whitelist)
   if(nrow(df)!=0){
      out[df$gene] <- 1
   }
   
   if(tag.features){
      names(out) <- paste0('gene_',mode,'.',names(out))
   }
   
   return(out)
}

# getGeneMutStatus(
#    driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv',
#    gene.whitelist.path=GENE_BIALLEL_WHITELIST,
#    mode='biallel'
# )
# 
# getGeneMutStatus(
#    driver.catalog.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180719_HMFregXXXXXXXX/sv-linx/XXXXXXXX.driver.catalog.tsv',
#    gene.whitelist.path=GENE_MONOALLEL_WHITELIST,
#    mode='monoallel'
# )











