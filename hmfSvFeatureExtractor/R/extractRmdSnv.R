#' Calculated regional mutation density (rmd)
#'
#' @param vcf.file Path to the vcf file
#' @param ref.genomeA  BSgenome reference genome. Default is
#' BSgenome.Hsapiens.UCSC.hg19. If another reference genome is indicated, it
#' will also need to be installed.
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
   vcf.file, ref.genome=DEFAULT_GENOME, verbose=F,
   sample.name=NULL, bin.size=1e6, as.matrix=T, ...
){
   #vcf.file='/Users/lnguyen//hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregCPCT_FR12246237_FR15412823_CPCT02380011/CPCT02380011T.purple.somatic.vcf.gz'

   df <- variantsFromVcf(vcf.file, ref.genome=ref.genome, verbose=verbose, ...)
   #df <- variantsFromVcf(vcf.file, ref.genome=ref.genome, verbose=verbose, vcf.filter='PASS')

   df_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(df)[1:4], df_colnames))){
      warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(df)[1:4] <- df_colnames
   }

   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]

   if(verbose){ message('Converting chrom name style to style in ref.genome...') }
   GenomeInfoDb::seqlevelsStyle(df$chrom)<- GenomeInfoDb::seqlevelsStyle(ref.genome)

   if(verbose){ message('Subsetting for SNVs...') }
   df <- df[nchar(df$ref)==1 & nchar(df$alt)==1,]

   ##----------------------------------------------------------------
   if(verbose){ message('Constructing genomic bins...') }
   chrom_lengths <- data.frame(
      chrom=BSgenome::seqinfo(ref.genome)@seqnames,
      length=BSgenome::seqinfo(ref.genome)@seqlengths
   )

   chrom_lengths <- chrom_lengths[ match(paste0('chr',1:22), chrom_lengths$chrom),]

   chrom_linear_offset <- (function(){
      v <- c(0, cumsum(as.numeric(chrom_lengths$length)))
      v <- v[-length(v)]
      names(v) <- chrom_lengths$chrom
      return(v)
   })()

   linearizeChromPos <- function(chrom, pos){
      #chrom=df$chrom
      #pos=df$pos
      chrom <- as.character(chrom)
      pos + unname(chrom_linear_offset[chrom])
   }

   genome_bins <- do.call(rbind, lapply(1:nrow(chrom_lengths), function(i){
      chrom_length <- chrom_lengths[i,'length']

      out <- data.frame( start=seq(1,chrom_length,by=bin.size) )
      out$end <- c(
         out$start[-1] - 1,
         chrom_length
      )

      out <- cbind(
         chrom=chrom_lengths[i,'chrom'],
         out
      )

      return(out)
   }))

   genome_bins$start_linear <- linearizeChromPos(genome_bins$chrom, genome_bins$start)
   genome_bins$end_linear <- linearizeChromPos(genome_bins$chrom, genome_bins$end)

   genome_bins$id <- 1:nrow(genome_bins)

   ##----------------------------------------------------------------
   if(verbose){ message('Counting muts per bin...') }
   counts <- structure(rep(0, nrow(genome_bins)), names=genome_bins$id )

   if(nrow(df)!=0){
      df$pos_linear <- linearizeChromPos(df$chrom, df$pos)

      df$bin <- unlist(lapply(df$pos_linear, function(i){
         #i=df$pos_linear[14509]
         if(is.na(i)){ return(NA) }
         genome_bins[
            genome_bins$start_linear <= i & i <= genome_bins$end_linear,
            'id'
         ]
      }))

      tab <- table(df$bin)
      counts[names(tab)] <- tab
   }

   if(verbose){ message('Returning output...') }


   if(!as.matrix){ return(unname(counts)) }

   out <- as.matrix(counts)
   rownames(out) <- paste0('r',rownames(out))

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
