#' Calculated regional mutation density (rmd)
#'
#' @param vcf.file Path to the vcf file
#' @param verbose Print progress messages?
#' @param sample.name If a character is provided, the header for the output
#' matrix will be named to this. If none is provided, the basename of the vcf
#' file will be used.
#' @param bin.size Size of genomic bins
#' @param as.matrix Return a matrix or vector?
#' @param ... Other arguments that can be passed to variantsFromVcf()
#'
#' @return A matrix or vector
#' @export
#'
extractRmdSnv <- function(
   vcf.file=NULL, df=NULL, sample.name=NULL, bin.size=1e6, as.matrix=T,
   verbose=F, ...
){
   #vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregCPCT_FR12246237_FR15412823_CPCT02380011/CPCT02380011T.purple.somatic.vcf.gz'
   #vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/200117_HMFregWIDE_FR17456804_FR17456800_WIDE01010504/WIDE01010504T.purple.somatic.vcf.gz'
   
   if(!is.null(vcf.file)){
      df <- vcfAsDataframe(vcf.file, verbose=verbose, ...)
   }
   #df <- vcfAsDataframe(vcf.file, verbose=verbose, vcf.filter='PASS')

   df_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(df)[1:4], df_colnames))){
      warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(df)[1:4] <- df_colnames
   }

   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]

   if(verbose){ message('Subsetting for SNVs...') }
   df <- df[nchar(df$ref)==1 & nchar(df$alt)==1,]
   
   if(nrow(df)!=0){
      if(verbose){ message('Converting chrom name style...') }
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- 'NCBI'
   }

   ##----------------------------------------------------------------
   if(verbose){ message('Counting muts per bin...') }
   genome_bins <- mkGenomeBins(bin.size=bin.size)
   GenomeInfoDb::seqlevelsStyle(genome_bins$chrom )<- 'NCBI'
   
   counts <- structure(rep(0, nrow(genome_bins)), names=genome_bins$bin_name )
   
   if(nrow(df)!=0){
      df <- df[df$chrom %in% genome_bins$chrom,]
      df$bin_name <- (function(){
         l <- Map(function(chrom,pos){
            which(
               chrom==genome_bins$chrom &
               pos>=genome_bins$start &
               pos<=genome_bins$end
            )
         }, df$chrom, df$pos)
         l[sapply(l,length)==0] <- NA
         genome_bins$bin_name[ unlist(l, use.names=F) ]
      })()
      
      tab <- table( df[!is.na(df$bin_name),'bin_name'] )
      counts[names(tab)] <- tab
   }
   
   ##----------------------------------------------------------------
   if(verbose){ message('Returning output...') }
   
   if(!as.matrix){ return(unname(counts)) }
   out <- as.matrix(counts)

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}

# m <- extractRmdSnv(
#    '/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregCPCT_FR12246237_FR15412823_CPCT02380011/CPCT02380011T.purple.somatic.vcf.gz',
#    vcf.filter='PASS', as.matrix=F, verbose=T
# )
