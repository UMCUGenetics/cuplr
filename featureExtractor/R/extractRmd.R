#' Calculated regional mutation density (rmd)
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
#' vcf.file
#' @param sample.name If a character is provided, the header for the output
#' matrix will be named to this. If none is provided, the basename of the vcf
#' file will be used.
#' @param genome.bins A dataframe containing predetermined bins. Must contain the columns: chrom, 
#' start, end, chrom_arm, start_linear, end_linear, bin_name
#' @param bin.size Alternatively if `genome.bins` is NULL, the user can specify the size of the
#' genomic bins
#' @param as.matrix Return a matrix or vector?
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param clonal.variants.only If TRUE, only variants with subclonal likelihood (from GRIDSS) < 0.8
#' will be selected
#' @param mut.types A character vector containing one or more of the following values: 'snv', 'dbs',
#' 'indel', or 'other'. Default is 'snv'. Only these mut types will be taken into account when
#' calculating the number of mutations per RMD bin
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to mutSigExtractor::variantsFromVcf()
#'
#' @return A matrix or vector
#' @export
#'
extractRmd <- function(
   vcf.file=NULL, df=NULL, sample.name=NULL, genome.bins=NULL, bin.size=1e6, as.matrix=T,
   vcf.filter='PASS', clonal.variants.only=T, mut.types='snv',
   verbose=F, ...
){
   # if(F){
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/200117_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/160825_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor//DO217817-from-jar//purple25/DO217817T.purple.somatic.vcf.gz'
   # 
   #    bin.size=1e6
   #    as.matrix=T
   #    vcf.filter='PASS'
   #    clonal.variants.only=T
   #    verbose=T
   # }
   
   if(!is.null(vcf.file)){
      df <- mutSigExtractor::variantsFromVcf(vcf.file, vcf.fields=c(1,2,4,5,7,8), vcf.filter=vcf.filter, verbose=verbose, ...)
   }
   
   if(clonal.variants.only & nrow(df)!=0){
      if(verbose){ message('Selecting clonal variants') }
      ## Split clonal/subclonal variants
      df$subclonal_prob <- as.numeric(mutSigExtractor::getInfoValues(df$info,'SUBCL'))
      df$subclonal_prob[is.na(df$subclonal_prob)] <- 0
      df$is_subclonal <- df$subclonal_prob >= 0.8
      df <- df[!df$is_subclonal,]
   }
   df$info <- NULL
   
   ## --------------------------------
   if(nrow(df)!=0){
      df_colnames <- c('chrom','pos','ref','alt')
      if(!(identical(colnames(df)[1:4], df_colnames))){
         warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
         colnames(df)[1:4] <- df_colnames
      }
      
      if(verbose){ message('Removing rows with multiple ALT sequences...') }
      df <- df[!grepl(',',df$alt),]
       
      if(verbose){ message('Converting chrom name style...') }
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- 'NCBI'
   }
   
   if(nrow(df)!=0){
      
      if(verbose){ message('Subsetting for: ',paste(mut.types, collapse=', '),'...') }
      
      df$mut_type <- 'other'
      df$mut_type[ nchar(df$ref)==1 & nchar(df$alt)==1 ] <- 'snv'
      df$mut_type[ nchar(df$ref)==2 & nchar(df$alt)==2 ] <- 'dbs'
      df$mut_type[ 
         ( nchar(df$ref)==1 & nchar(df$alt)>1 ) |
         ( nchar(df$ref)>1 & nchar(df$alt)==1 )
      ] <- 'indel'
      
      df <- df[df$mut_type %in% mut.types,]
   }

   ## --------------------------------
   if(verbose){ message('Counting muts per bin...') }
   if(is.null(genome.bins)){
      genome.bins <- mkGenomeBins(bin.size=bin.size)
   }
   GenomeInfoDb::seqlevelsStyle(genome.bins$chrom )<- 'NCBI'
   
   counts <- structure(rep(0, nrow(genome.bins)), names=genome.bins$bin_name )
   
   if(nrow(df)!=0){
      df <- df[df$chrom %in% genome.bins$chrom,]
      df$bin_name <- (function(){
         #counter <- 0
         l <- Map(function(chrom,pos){
            #counter <<- counter + 1
            #print(counter)
            which(
               chrom==genome.bins$chrom &
               pos>=genome.bins$start &
               pos<=genome.bins$end
            )
         }, df$chrom, df$pos)
         
         item_lengths <- sapply(l,length)
         
         if(any(item_lengths>1)){ 
            stop("Some variants fall into multiple bins (overlapping bin coordinates)") 
         }
         
         l[item_lengths==0] <- NA
         genome.bins$bin_name[ unlist(l, use.names=F) ]
      })()
      
      tab <- table( df[!is.na(df$bin_name),'bin_name'] )
      counts[names(tab)] <- tab
   }
   
   ## --------------------------------
   if(verbose){ message('Returning output...') }
   
   if(!as.matrix){ return(counts) }
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
