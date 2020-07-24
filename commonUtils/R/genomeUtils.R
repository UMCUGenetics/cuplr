## hg19 constants
## centromere positions
CENTRO_POS_HG19 <- c(
   chr1=123035434,chr2=93826171,chr3=92004854,chr4=51160117,chr5=47905641,
   chr6=60330166,chr7=59554331,chr8=45338887,chr9=48867679,chr10=40754935,
   chr11=53144205,chr12=36356694,chr13=17500000,chr14=17500000,chr15=18500000,
   chr16=36835801,chr17=23763006,chr18=16960898,chr19=26181782,chr20=27869569,
   chr21=12788129,chr22=14500000,chrX=60132012,chrY=11604553
)

CHROM_LENGTHS_HG19 <- c(
   chr1=249250621,chr2=243199373,chr3=198022430,chr4=191154276,chr5=180915260,
   chr6=171115067,chr7=159138663,chr8=146364022,chr9=141213431,chr10=135534747,
   chr11=135006516,chr12=133851895,chr13=115169878,chr14=107349540,chr15=102531392,
   chr16=90354753,chr17=81195210,chr18=78077248,chr19=59128983,chr20=63025520,
   chr21=48129895,chr22=51304566,chrX=155270560,chrY=59373566,chrM=16571
)

####################################################################################################
#' Find overlapping intervals
#'
#' @param start1 A vector of start positions (#1)
#' @param end1 A vector of end positions (#1)
#' @param start2 A vector of start positions (#2)
#' @param end2 A vector of end positions (#2)
#' @param verbose Show progress?
#'
#' @return A list of indexes of #1 for which #2 overlaps
#' @export
#'
isOverlapping <- function(start1, end1, start2, end2, verbose=F){
   # start1=genome_bins$start_linear
   # end1=genome_bins$end_linear
   #
   # start2=linearizeChromPos(cnv$chromosome, cnv$start)
   # end2=linearizeChromPos(cnv$chromosome, cnv$end)

   if( length(start1)!=length(end1) ){
      stop('start1 and end1 must be the same length')
   }
   if( length(start2)!=length(end2) ){
      stop('start2 and end2 must be the same length')
   }

   if(verbose){ pb <- txtProgressBar(max=length(start1), style=3) }
   lapply(1:length(start1), function(i){
      #i=1
      # out <- which(pmax(start1[i],start2) <= pmin(end1[i], end2))
      # if(length(out)!=0){ out } else { NA }
      if(verbose){ setTxtProgressBar(pb, i) }
      which(pmax(start1[i],start2) <= pmin(end1[i], end2))
   })
}

#' Find overlapping genomic intervals
#'
#' @param chrom1 A vector or chromosome names (#1)
#' @param start1 A vector of start positions (#1)
#' @param end1 A vector of end positions (#1)
#' @param chrom2 A vector or chromosome names (#1)
#' @param start2 A vector of start positions (#2)
#' @param end2 A vector of end positions (#2)
#'
#' @return A list of indexes of #1 for which #2 overlaps
#' @export
#'
isOverlappingChromPos <- function(
   chrom1=NULL, start1=NULL, end1=NULL, chrom2=NULL, start2=NULL, end2=NULL,
   df1=NULL, df2=NULL
){
   ## Debug
   # chrom1=bed$chrom
   # start1=bed$start
   # end1=bed$end
   #
   # chrom2=cnv$chrom
   # start2=cnv$start
   # end2=cnv$end
   #
   # df1 <- data.frame(chrom=chrom1, start=start1, end=end1, stringsAsFactors=F)
   # df2 <- data.frame(chrom=chrom2, start=start2, end=end2, stringsAsFactors=F)
   # colnames(df1)[1:3] <- colnames(df2)[1:3] <- c('chrom','start','end')

   if(!is.null(df1)){
      chrom1 <- df1[,1]
      start1 <- df1[,2]
      end1 <- df1[,3]
   }

   if(!is.null(df2)){
      chrom2 <- df2[,1]
      start2 <- df2[,2]
      end2 <- df2[,3]
   }

   GenomeInfoDb::seqlevelsStyle(chrom1)<- 'NCBI'
   GenomeInfoDb::seqlevelsStyle(chrom2)<- 'NCBI'

   chrom_lookup <- as.factor(unique(c(chrom1,chrom2)))

   chrom1 <- factor(chrom1, levels=chrom_lookup)
   chrom2 <- factor(chrom2, levels=chrom_lookup)

   pad_width <- nchar(max(c(start1, end1, start2, end2)))

   chromPosAsNumeric <- function(chrom, pos){
      as.numeric(paste0(
         as.integer(chrom),
         formatC(pos, width=pad_width, format='d', flag='0')
      ))
   }

   start1 <- chromPosAsNumeric(chrom1, start1)
   end1 <- chromPosAsNumeric(chrom1, end1)

   start2 <- chromPosAsNumeric(chrom2, start2)
   end2 <- chromPosAsNumeric(chrom2, end2)

   isOverlapping(start1, end1, start2, end2)
}

####################################################################################################
#' Get chromosome arm
#'
#' @param df A dataframe containing the columns: chrom, start, end
#' @param centro.pos A named integer vector of the centromere positions (default: hg19)
#' @param chrom.arm.names A character vector in the form c('1p','1q','2p','2q', ...). The default
#' 'auto' means that the human chromosome arm names are used. Note that chroms 13, 14, 15, 21, 22
#' are considered to only have the long (i.e. q) arm.
#' @param arm.only Only return the arm and not the chrom number?
#' @param show.warnings Show warning messages?
#'
#' @return A character vector of chromosome arm names
#' @export
#'
getChromArm <- function(
   df=NULL, centro.pos=CENTRO_POS_HG19, one.armed.chroms=c(13,14,15,21,22),
   arm.only=F, seq.levels.style='NCBI', show.warnings=F
){

   colnames(df)[1:3] <- c('chrom','start','end')

   df$chrom <- as.character(df$chrom)
   GenomeInfoDb::seqlevelsStyle(df$chrom)<- GenomeInfoDb::seqlevelsStyle(names(centro.pos))[1]

   one.armed.chroms <- as.character(one.armed.chroms)
   GenomeInfoDb::seqlevelsStyle(one.armed.chroms)<- GenomeInfoDb::seqlevelsStyle(names(centro.pos))[1]

   df$centro_pos <- centro.pos[ df$chrom ]

   df$is_pq <- df$start < df$centro_pos & df$end > df$centro_pos
   df$is_p <- df$start < df$centro_pos & df$end < df$centro_pos
   df$is_q <- df$start >= df$centro_pos & df$end >= df$centro_pos

   if(any(df$is_pq) & show.warnings){
      warning(sum(df$is_pq)," intervals overlap the centromere. Assigning as 'pq'")
   }
   df$arm <- c('pq','p','q')[
      max.col(df[,c('is_pq','is_p','is_q')])
   ]

   df$arm[df$chrom %in% one.armed.chroms] <- 'q'

   if(!arm.only){
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- seq.levels.style
      paste0(df$chrom,df$arm)
   } else {
      df$arm
   }
}

####################################################################################################
#' Get cytoband
#'
#' @param df A dataframe containing the columns: chrom, start, end
#' @param cytobands.path Path to the cytobands txt file with the columns: chrom, start, end, band, stain
#' @param greedy.multi.match If TRUE, when multiple cytobands match to a region, only the first
#' cytoband will be reported
#' @param verbose Show progress?
#'
#' @return A character vector
#' @export
#'
getCytoband <- function(df, cytobands.path=CYTOBANDS, greedy.multi.match=F, verbose=F){
   #df=read.delim('/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/160704_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv')
   #df=df[,1:3]
   colnames(df)[1:3] <- c('chrom','start','end')

   cytobands <- read.delim(cytobands.path)
   cytobands$start <- cytobands$start + 1 ## Make start 1-based

   ## Calculate linear coords
   df$start_linear <- linearizeChromPos(df$chrom, df$start)
   df$end_linear <- linearizeChromPos(df$chrom, df$end)
   cytobands$start_linear <- linearizeChromPos(cytobands$chrom, cytobands$start)
   cytobands$end_linear <- linearizeChromPos(cytobands$chrom, cytobands$end)

   ##
   overlaps <- isOverlapping(
      df$start_linear, df$end_linear,
      cytobands$start_linear, cytobands$end_linear,
      verbose=verbose
   )

   if(!greedy.multi.match){
      unlist(lapply(overlaps, function(i){
         if(length(i)==0){ return(NA) }
         paste( cytobands[i,'band'], collapse=';')
      }), use.names=F)
   } else {
      unlist(lapply(overlaps,function(i){
         if(length(i)==0){ return(NA) }
         cytobands[i[1],'band']
      }))

   }
}

####################################################################################################
#' Converts chrom/pos into linear coordinates
#'
#' @param chrom A character vector indicating the chromosome
#' @param pos An integer vector indicating the position
#' @param df Alternative input to chrom and pos. A dataframe with chrom in 1st col and pos in 2nd col
#' @param chrom.lengths A named vector of chromosome lengths (default: hg19)
#'
#' @return An integer vector of linear coordinates
#' @export
#'
linearizeChromPos <- function(chrom=NULL, pos=NULL, df=NULL, chrom.lengths=CHROM_LENGTHS_HG19){

   if(!is.null(df)){
      chrom <- df[,1]
      pos <- df[,2]
   }

   if(length(chrom)!=length(pos)){ stop('`chrom` and `pos` are not the same length') }

   GenomeInfoDb::seqlevelsStyle(chrom)<- GenomeInfoDb::seqlevelsStyle(names(chrom.lengths))[1]

   offsets <- cumsum(as.numeric(chrom.lengths))
   offsets <- c(0,offsets[-length(offsets)])
   names(offsets) <- names(chrom.lengths)

   unname(pos + offsets[chrom])
}

#' Make a look up table of genomic bin intervals
#'
#' @param bin.size Size of each genome bin
#' @param keep.chroms A character vector specifying which chromosomes to keep (chromosome names
#' should be in the style of the vcf). To keep autosomal and sex chromosomes for example use:
#' keep.chroms=c(1:22,'X','Y')
#' @param split.centro.bins When bins overlap centromere positions, split these into 2 bins?
#' @param centro.pos A named integer vector of the centromere positions (default: hg19)
#' @param chrom.lengths A named vector of chromosome lengths (default: hg19)
#'
#' @return A dataframe containing the intervals of each bin (in chr:start:end and
#' linear_start:linear_end form)
#' @export
#'
mkGenomeBins <- function(
   bin.size=1E6, keep.chroms=c(1:22,'X'), split.centro.bins=T,
   centro.pos=CENTRO_POS_HG19, chrom.lengths=CHROM_LENGTHS_HG19
){

   keep.chroms <- as.character(keep.chroms)
   GenomeInfoDb::seqlevelsStyle(keep.chroms)<- GenomeInfoDb::seqlevelsStyle(names(chrom.lengths))[1]
   GenomeInfoDb::seqlevelsStyle(names(centro.pos))<- GenomeInfoDb::seqlevelsStyle(names(chrom.lengths))[1]

   chrom.lengths <- chrom.lengths[ match(keep.chroms, chrom.lengths$chrom),]

   genome_bins <- do.call(rbind, lapply(names(chrom.lengths), function(i){
      chrom_length <- chrom.lengths[i]

      out <- data.frame( start=seq(1,chrom_length,by=bin.size) )
      out$end <- c(
         out$start[-1] - 1,
         chrom_length
      )

      out <- cbind(chrom=i, out)

      return(out)
   }))

   if(split.centro.bins){
      genome_bins$chrom_arm <- getChromArm(
         genome_bins[,c('chrom','start','end')],
         show.warnings=F
      )

      pq_index <- grep('pq$',genome_bins$chrom_arm)
      for(i in 1:length(pq_index)){
         #i=pq_pos[1]
         current_index <- pq_index[i]
         genome_bins <- insertRow(genome_bins, row.num=current_index)
         centro_pos <- centro.pos[ genome_bins[current_index,'chrom'] ]

         genome_bins[current_index,'end'] <- centro_pos-1
         genome_bins[current_index+1,'start'] <- centro_pos

         pq_index <- pq_index+1
      }

      genome_bins$chrom_arm <- getChromArm(
         genome_bins[,c('chrom','start','end')],
         show.warnings=F
      )

      #genome_bins[grep('pq',genome_bins$chrom_arm),]
      #genome_bins$chrom_arm <- NULL
   }

   genome_bins$start_linear <- linearizeChromPos(df=genome_bins[,c('chrom','start')])
   genome_bins$end_linear <- linearizeChromPos(df=genome_bins[,c('chrom','end')])

   genome_bins$id <- 1:nrow(genome_bins)

   return(genome_bins)
}




