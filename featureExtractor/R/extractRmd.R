#' Calculated regional mutation density (rmd)
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
#' vcf.file
#' @param sample.name If a character is provided, the header for the output
#' matrix will be named to this. If none is provided, the basename of the vcf
#' file will be used.
#' @param genome.bins A dataframe containing predetermined bins. Must contain the columns: chrom, 
#' start, end, bin_name. Alternatively, if a path to a txt (as a string) is provided, the txt will
#' be loaded.
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
   vcf.file=NULL, df=NULL, sample.name=NULL, genome.bins=RMD_BINS_PATH, as.matrix=T,
   vcf.filter='PASS', clonal.variants.only=F, mut.types='snv',
   verbose=F, ...
){
   # if(F){
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/200117_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/160825_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   #    vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor//DO217817-from-jar//purple25/DO217817T.purple.somatic.vcf.gz'
   #    
   #    bin.size=1e6
   #    genome.bins=rmd_bins.hg38
   #    genome.bins=rmd_bins.hg19
   #    as.matrix=T
   #    vcf.filter='PASS'
   #    clonal.variants.only=F
   #    mut.types='snv'
   #    verbose=T
   # }
   
   ## Init --------------------------------
   if(!is.null(vcf.file)){
      df <- mutSigExtractor::variantsFromVcf(vcf.file, vcf.fields=c(1,2,4,5,7,8), vcf.filter=vcf.filter, verbose=verbose, ...)
      #df <- mutSigExtractor::variantsFromVcf(vcf.file, vcf.fields=c(1,2,4,5,7,8), vcf.filter=vcf.filter, verbose=verbose)
   }
   
   ## Subset for clonal variants
   if(clonal.variants.only & nrow(df)!=0){
      if(verbose){ message('Selecting clonal variants') }
      ## Split clonal/subclonal variants
      df$subclonal_prob <- as.numeric(mutSigExtractor::getInfoValues(df$info,'SUBCL'))
      df$subclonal_prob[is.na(df$subclonal_prob)] <- 0
      df$is_subclonal <- df$subclonal_prob >= 0.8
      df <- df[!df$is_subclonal,]
   }
   df$info <- NULL
   
   ## Sanitize vcf input
   if(nrow(df)!=0){
      df_colnames <- c('chrom','pos','ref','alt')
      if(!(identical(colnames(df)[1:4], df_colnames))){
         warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
         colnames(df)[1:4] <- df_colnames
      }
      
      if(verbose){ message('Removing rows with multiple ALT sequences...') }
      df <- df[!grepl(',',df$alt),]
      
      if(verbose){ message("Setting chrom name style for vcf to 'NCBI'...") }
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- 'NCBI'
   }
   
   ## Subset for mut types
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
   
   ## Initialize counts
   if(is.character(genome.bins)){
      if(verbose){ message("Reading in genome bins @ ", genome.bins) }
      genome.bins <- read.delim(genome.bins) 
   }
   
   counts <- structure(
      rep(0, nrow(genome.bins)), 
      names=naturalsort::naturalsort(unique(genome.bins$bin_name))
   )
   
   ## Main
   if(nrow(df)!=0){
      
      ##
      if(verbose){ message("Setting chrom name style for `genome.bins` to 'NCBI'...") }
      genome.bins$chrom <- as.character(genome.bins$chrom)
      GenomeInfoDb::seqlevelsStyle(genome.bins$chrom)<- 'NCBI'
      
      ##
      if(verbose){ message("Convert chrom:pos to integer representation...") }
      df <- df[df$chrom %in% genome.bins$chrom,]
      uniq_chroms <- naturalsort::naturalsort(unique(genome.bins$chrom))
      genome.bins$chrom <- factor(genome.bins$chrom, uniq_chroms)
      df$chrom <- factor(df$chrom, uniq_chroms)
      
      max_pos_digits <- nchar(max(c(genome.bins$start, genome.bins$end, df$pos)))
      chromPosToInt <- function(chrom, pos, n.digits=max_pos_digits){
         # chrom=genome.bins$chrom
         # pos=genome.bins$start
         # n.digits=max_pos_digits
         
         pos_char <- formatC(pos, width=n.digits, format='d', flag='0')
         chrom_int <- as.integer(chrom)
         as.numeric(paste0(chrom_int, pos_char))
      }
      
      genome.bins$start_int <- chromPosToInt(genome.bins$chrom, genome.bins$start)
      genome.bins$end_int <- chromPosToInt(genome.bins$chrom, genome.bins$end)
      df$pos_int <- chromPosToInt(df$chrom, df$pos)
      
      ##
      if(verbose){ message("Determining RMD bin for each mutation in the vcf...") }
      matched_bins <- lapply(df$pos_int, function(i){
         #i=df$pos_int[8221]
         genome.bins$bin_name[ i>=genome.bins$start_int & i<=genome.bins$end_int ]
      })
      
      match_lengths <- sapply(matched_bins, length)
      n_warn_matches <- sum(match_lengths>1 | match_lengths==0)
      
      if(n_warn_matches>0){ 
         warn_matches_summary <- table(match_lengths)
         warn_matches_summary <- warn_matches_summary[names(warn_matches_summary)!=1]
         warn_matches_summary <- paste0(warn_matches_summary, ' variants in ', names(warn_matches_summary),' bins')
         warn_matches_summary <- paste(warn_matches_summary, collapse=', ')
         
         warning(
            n_warn_matches,'/',nrow(df)," variants were matched to zero or multiple bins:\n",
            warn_matches_summary,'.\n',
            "For variants matched to >=2 bins, the first bin will be assigned"
         ) 
      }
      
      if(verbose){ message("Counting bin occurrences...") }
      matched_bins[match_lengths==0] <- NA
      df$bin_name <- sapply(matched_bins,`[[`,1)
      tab <- table( df$bin_name[!is.na(df$bin_name)] )
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

