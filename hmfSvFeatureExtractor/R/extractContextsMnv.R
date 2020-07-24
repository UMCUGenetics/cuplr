#' Extract MNV contexts
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to 
#' vcf.file
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param keep.chroms A character vector specifying which chromosomes to keep (chromosome names
#' should be in the style of the vcf). To keep autosomal and sex chromosomes for example use:
#' keep.chroms=c(1:22,'X','Y')
#' @param len.cap MNVs longer than this length will be put into the same bin. MNV length is defined 
#' as max(nchar(REF),nchar(ALT))
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param tag.features If TRUE, will tag the value names with the feature type
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to mutSigExtractor::variantsFromVcf()
#'
#' @return An integer vector
#' @export
#'
extractContextsMnv <- function(
   vcf.file=NULL, df=NULL, vcf.filter='PASS', keep.chroms=c(1:22,'X'),
   len.cap=5, ref.genome=mutSigExtractor::DEFAULT_GENOME,
   tag.features=T, verbose=F, ...
){
   #vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/180424_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/160709_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/190723_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   
   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- mutSigExtractor::variantsFromVcf(
         vcf.file, vcf.filter=vcf.filter, keep.chroms=keep.chroms,
         ref.genome=ref.genome, verbose=verbose,
         ...
      )
      
      # df <- mutSigExtractor::variantsFromVcf(
      #    vcf.file, vcf.filter=vcf.filter,
      #    ref.genome=ref.genome, verbose=verbose
      # )
   }
   
   if(verbose){ message('Counting MNV contexts...') }
   context_names <- c(
      paste0('mnv_del.len.',3:len.cap),
      paste0('mnv_ins.len.',3:len.cap),
      paste0('mnv_neutral.len.',3:len.cap)
   )
   
   counts <- structure(rep(0, length(context_names)), names=context_names)
   
   if(is.data.frame(df)){
      df <- df[!grepl(',',df$alt),]
      
      df$ref_len <- nchar(df$ref)
      df$alt_len <- nchar(df$alt)
      
      ## Subset for MNVs (excluding DBSs)
      df <- df[
         with(df,{ (ref_len>=3 & alt_len>=2) | (ref_len>=2 & alt_len>=3) })
      ,]
      
      df$max_len <- pmax(df$ref_len, df$alt_len)
      df$max_len[df$max_len>len.cap] <- len.cap
      
      df$type <- 'neutral'
      df$type[ df$ref_len > df$alt_len ] <- 'del'
      df$type[ df$ref_len < df$alt_len ] <- 'ins'
      
      df$context <- with(df,{
         paste0('mnv_',type,'.len.',max_len)
      })
      
      tab <- table(df$context)
      counts[names(tab)] <- tab
   }
   
   if(tag.features){
      names(counts) <- paste0('mnv.',names(counts))
   }
   
   return(counts)
}

