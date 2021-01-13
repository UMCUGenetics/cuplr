#' Extract SNV/indel contexts/signatures
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to 
#' vcf.file
#' @param sample.name If a character is provided, the header for the output
#' matrix will be named to this. If none is provided, the basename of the vcf
#' file will be used.
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param keep.chroms A character vector specifying which chromosomes to keep (chromosome names
#' should be in the style of the vcf). To keep autosomal and sex chromosomes for example use:
#' keep.chroms=c(1:22,'X','Y')
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param clonal.variants.only If TRUE, only variants with subclonal likelihood (from GRIDSS) < 0.8
#' will be selected
#' @param as.matrix Return a matrix or vector?
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to mutSigExtractor::variantsFromVcf()
#'
#' @return An matrix or integer vector
#' @export
#'
extractContextsSmnvIndel <- function(
   vcf.file=NULL, df=NULL, sample.name=NULL, vcf.filter='PASS', keep.chroms=c(1:22,'X'),
   ref.genome=mutSigExtractor::DEFAULT_GENOME, clonal.variants.only=T,
   as.matrix=T, verbose=F, ...
){

   if(F){
      vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
      
      vcf.filter='PASS'
      keep.chroms=c(1:22,'X')
      
      ref.genome=mutSigExtractor::DEFAULT_GENOME
      tag.features=T
      verbose=T
      
      clonal.variants.only=T
   }

   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- mutSigExtractor::variantsFromVcf(
         vcf.file, vcf.fields=c(1,2,4,5,7,8),
         vcf.filter=vcf.filter, keep.chroms=keep.chroms,
         ref.genome=ref.genome, verbose=verbose,
         ...
      )
      
      # df <- mutSigExtractor::variantsFromVcf(
      #    vcf.file, vcf.fields=c(1,2,4,5,7,8),
      #    vcf.filter=vcf.filter, keep.chroms=keep.chroms,
      #    ref.genome=ref.genome, verbose=verbose
      # )
      
      if(clonal.variants.only){
         if(verbose){ message('Selecting clonal variants') }
         ## Split clonal/subclonal variants
         df$subclonal_prob <- as.numeric(getInfoValues(df$info,'SUBCL')[,1])
         df$subclonal_prob[is.na(df$subclonal_prob)] <- 0
         df$is_subclonal <- df$subclonal_prob >= 0.8
         df <- df[!df$is_subclonal,]
      }
      df$info <- NULL
   }
   
   df_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(df)[1:4], df_colnames))){
      warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(df)[1:4] <- df_colnames
   }
   
   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]
   
   if(nrow(df)!=0){
      if(verbose){ message('Converting chrom name style...') }
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- 'NCBI'
   }
   
   ##----------------------------------------------------------------
   l <- list()
   
   if(verbose){ message('\n## Extracting SNV contexts...') }
   l$snv <- mutSigExtractor::extractSigsSnv(df=df, output='contexts',verbose=verbose)[,1]
   
   if(verbose){ message('\n## Extracting indel contexts...') }
   l$indel <- mutSigExtractor::extractSigsIndel(df=df, output='contexts', method='PCAWG', verbose=verbose)[,1]
   
   if(verbose){ message('\n## Extracting DBS contexts...') }
   l$dbs <- mutSigExtractor::extractSigsDbs(df=df, output='contexts', verbose=verbose)[,1]
   
   if(verbose){ message('\n## Counting MNVs...') }
   ## excluding DBSs
   df$ref_len <- nchar(df$ref)
   df$alt_len <- nchar(df$alt)
   l$mnv <- c(
      counts=sum(with(df,{ (ref_len>=3 & alt_len>=2) | (ref_len>=2 & alt_len>=3) }))
   )
   
   ##----------------------------------------------------------------
   if(verbose){ message('Returning output...') }
   out <- do.call(c, l)
   
   if(!as.matrix){ return(unname(out)) }
   out <- as.matrix(out)
   
   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }
   
   return(out)
}

