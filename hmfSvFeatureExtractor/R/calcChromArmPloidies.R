#' Calculate overall chrom arm copy numbers
#'
#' @description This function first rounds copy numbers (CN) to integers so that CN segments can be
#'   grouped together. Per chrom arm, the coverage of each CN category is calculated (i.e.
#'   cumulative segment size). The chrom arm CN is (roughly) defined as the CN category with the
#'   highest cumulative segment size
#'
#' @param cnv.file Path to purple cnv file
#' @param out.file Path to output file. If NULL, returns a named vector
#' @param min.rel.cum.segment.size If a chrom arm has a CN category that covers >0.5 (i.e 50%;
#'   default) of a chrom arm, this CN is the copy number of the arm
#' @param max.rel.cum.segment.size.diff This value (default 0.1) determines whether which CN
#'   categories are considered to cover equal lengths of the chrom arm. For example, (by default) 2
#'   CN categories covering 0.40 and 0.31 of a chrom arm are considered equally contributing. When
#'   these CN categories have similar cumulative segment size as the one with the highest, if one of
#'   these have the same CN as the genome CN, return the genome CN. Otherwise, simply return the one
#'   with the highest segment support (as is done above).
#' @param chrom.arm.names A character vector in the form c('1p','1q','2p','2q', ...). The default
#'   'auto' means that the human chromosome arm names are used. Note that chroms 13, 14, 15, 21, 22
#'   are considered to only have the long (i.e. q) arm.
#' @param chrom.arm.split.method Which method to determine the chromosome arm coords? If 'hmf', uses
#'   'method' column from purple cnv file to determine centromere positions (i.e. p/q arm split
#'   point). If 'universal', uses the a (processed) gap.txt.gz table from the UCSC genome browser to
#'   determine centromere positions. These 2 methods should in theory be identical, unless the HMF
#'   pipeline code changes.
#' @param tag.features If TRUE, will tag the value names with the feature type
#' @param verbose Show progress messages?
#'
#' @return A named vector of chrom arm copy numbers, or specified writes a table to out.file if
#'   specified
#' @export
#'
#' @examples
calcChromArmPloidies <- function(
   cnv.file, out.file=NULL, 
   min.rel.cum.segment.size=0.5, max.rel.cum.segment.size.diff=0.1,
   chrom.arm.split.method='universal', 
   centromere.positions.path=CENTROMERE_POSITIONS, 
   one.armed.chroms=c(13,14,15,21,22),
   chrom.arm.names='auto', ignore.chroms=NULL,
   tag.features=T, verbose=T
){
   #cnv.file="/Users/lnguyen//hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/PCAWG_2020/vcf/somatic/cnv/cna_annotated/005e85a3-3571-462d-8dc9-2babfc7ace21.consensus.20170119.somatic.cna.annotated.txt"
   
   cnv <- read.delim(cnv.file, check.names=F, stringsAsFactors=F)
   colnames(cnv) <- sub('#','',colnames(cnv))
   
   if(!is.null(ignore.chroms)){
      cnv <- cnv[!(cnv$chrom %in% ignore.chroms),]
   }
   
   #--------- Pre-calculations ---------#
   # if(chrom.arm.split.method=='hmf'){
   #    cnv <- cnv[,c('chromosome','start','end','copyNumber','segmentStartSupport','segmentEndSupport','method')]
   # }
   colnames(cnv)[1:4] <- c('chrom','start','end','copy_number')
   cnv$chrom <- gsub('chr','',cnv$chrom)
   
   cnv$segment_size <- (cnv$end - cnv$start) + 1
   cnv$copy_number[cnv$copy_number < 0] <- 0 ## Ceiling negative values
   cnv$copy_number_int <- round(cnv$copy_number)
   
   #--------- Split by chromosome and arm ---------#
   if(verbose){ message('Splitting CN segments by chromosome...') }
   cnv_split_pre <- split(cnv, cnv$chrom)
   cnv_split_pre <- cnv_split_pre[unique(cnv$chrom)] ## maintain chrom order
   
   if(verbose){ message('Splitting CN segments by arm...') }
   if(chrom.arm.split.method=='hmf'){
      cnv_split <- lapply(cnv_split_pre,function(i){
         #i <- cnv_split_pre[[22]]
         if(!('LONG_ARM' %in% i$method)){ ## One armed chromosomes are marked by LONG_ARM in method
            p_arm_end_index <- which(i$segmentEndSupport=='CENTROMERE')
            l <- list(
               p=i[1:p_arm_end_index,],
               q=i[(p_arm_end_index+1):nrow(i),]
            )
         } else {
            ## Consider one armed chroms as only having q arm
            l <- list(
               q=i
            )
         }
         names(l) <- paste0(i[1,'chrom'],names(l)) ## Make list names in the form 1p,1q,etc before flattening list
         return(l)
      })
   }
   
   if(chrom.arm.split.method=='universal'){
      centro_pos <- read.delim(centromere.positions.path, stringsAsFactors=F)
      centro_pos <- structure(centro_pos$pos, names=centro_pos$chrom)
      
      cnv_split <- lapply(cnv_split_pre,function(i){
         #i <- cnv_split_pre[['11']]
         chrom <- i[1,'chrom']
         if(!(chrom %in% one.armed.chroms)){
            q_arm_start_index <- which.min(i$end < centro_pos[chrom])
            
            l <- list(
               p=i[1:(q_arm_start_index-1),],
               q=i[q_arm_start_index:nrow(i),]
            )
         } else {
            l <- list(
               q=i
            )
         }
         names(l) <- paste0(i[1,'chrom'],names(l)) ## Make list names in the form 1p,1q,etc before flattening list
         return(l)
      })
   }
   
   ## Flatten 2-level list into 1-level list
   names(cnv_split) <- NULL
   cnv_split <- do.call(c, cnv_split)
   rm(cnv_split_pre)
   
   #--------- Calc arm ploidies ---------#
   if(verbose){ message('Calculating preliminary arm ploidies...') }
   cn_segment_support <- lapply(cnv_split, function(i){
      #i=cnv_split$`Yp`
      df <- aggregate(i$segment_size, by=list(i$copy_number_int), FUN=sum)
      colnames(df) <- c('copy_number_int','cum_segment_size')
      df <- as.data.frame(lapply(df, as.integer)) ## aggregate can return lists instead of vectors as output
      
      df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)
      
      #df[which.max(df$cum_segment_size_rel),]
      if(nrow(df)!=0){
         df <- df[order(df$cum_segment_size_rel, decreasing=T),]
         return(df)
      } else {
         df[1,] <- NA
         return(df)
      }
      
   })
   
   ## CN with most frequent total segment support is preliminary CN
   arm_cn_prelim <- unlist(lapply(cn_segment_support, function(i){ i[1,'copy_number_int'] }))
   
   if( sum(!is.na(arm_cn_prelim)) / length(arm_cn_prelim) < 0.6 ){
      if(verbose){ message('Warning: CN data contains too many NAs. Returning a vector of NAs') }
      genome_cn <- NA ## Genome ploidy cannot be calculated if there is too little CN data
   } else {
      genome_cn <- as.numeric(names(
         sort(table(arm_cn_prelim),decreasing=T)
      )[1])
   }
   
   if(is.na(genome_cn)){
      arm_cn <- structure(rep(NA, length(cn_segment_support)), names=names(cn_segment_support))
   } else {
      if(verbose){ message('Calculating final arm ploidies...') }
      arm_cn <- unlist(lapply(cn_segment_support, function(i){
         #i=cn_segment_support[['12p']]
         
         if(is.na(i[1,1])){ return(NA) }
         
         ## E.g. >=50% of arm has CN of 2, then this is the CN
         if( i[1,'cum_segment_size_rel'] >= min.rel.cum.segment.size ){
            return(i[1,'copy_number_int'])
         }
         
         #' When multiple CNs have similar segment support as the one with the highest, if one 
         #' of these have the same CN as the genome CN, return the genome CN. Otherwise, simply return 
         #' the one with the highest segment support (as is done above)
         i$diffs <- i[1,'cum_segment_size_rel'] - i[,'cum_segment_size_rel']
         cn_doubt <- i[i$diffs < max.rel.cum.segment.size.diff,'copy_number_int']
         
         if(any(cn_doubt==genome_cn)){
            return(genome_cn)
         } else {
            return(i[1,'copy_number_int'])
         }
      }))
   }
   
   ploidies <- c(arm_cn, genome=genome_cn)
   
   #--------- Ensure consistent output (e.g. for when CN data is missing) ---------#
   if(verbose){ message('Preparing output...') }
   if(is.null(chrom.arm.names)){
      return(ploidies)
   }
   
   if(chrom.arm.names=='auto'){
      chrom_arm_names <- paste0(rep(c(1:22,'X'),each=2),c('p','q'))
      chrom_arm_names <- chrom_arm_names[!(chrom_arm_names %in% paste0(one.armed.chroms,'p'))] ## Keep only q arm for one arm chromosomes
      chrom_arm_names <- c(chrom_arm_names,'genome')
   } else {
      chrom_arm_names <- c(chrom.arm.names,'genome')
   }
   
   out <- structure(rep(NA,length(chrom_arm_names)), names=chrom_arm_names)
   ploidies <- ploidies[names(ploidies) %in% chrom_arm_names] ## Rm unwanted chroms/arms (e.g. mitochondria)
   out[names(ploidies)] <- ploidies
   
   if(tag.features){
      names(out) <- paste0('ploidy.',names(out))
   }
   
   if(is.null(out.file)){ 
      return(out) 
   } else {
      out <- data.frame(chrom=names(out),ploidy=out,row.names=NULL)
      if(verbose){ message('Writing output...') }
      write.tsv(out, out.file)
   }
   
}

################################################################################
#' Calculate change in chrom arm copy number compared to genome copy number
#'
#' @param x Output from calcChromArmPloidies() as a matrix or vector
#'
#' @return A matrix or vector
#' @export
#'
calcChromArmCnChange <- function(x){
   if(is.matrix(x) | is.data.frame(x)){
      arm <- x[,-ncol(x)]
      genome <- x[,ncol(x)]
   } else if(is.numeric(x)){
      arm <- x[-length(x)]
      genome <- x[length(x)]
   } else {
      stop('`x` must be a vector or matrix')
   }
   
   main <- function(x, mode){
      x_diff <- arm - genome
      if(mode=='gain'){ x_diff[x_diff < 0] <- 0 } 
      else { x_diff[x_diff > 0] <- 0 }
      x_diff <- abs(x_diff)
      
      if(is.matrix(x) | is.data.frame(x)){
         colnames(x_diff) <- gsub('^\\w+[.]',paste0(mode,'.'), colnames(x_diff))
      } else {
         names(x_diff) <- gsub('^\\w+[.]',paste0(mode,'.'), names(x_diff))
      }
      
      return(x_diff)
   }
   
   if(is.matrix(x) | is.data.frame(x)){
      cbind(
         main(x, 'gain'),
         main(x, 'loss')
      )
   } else {
      c(
         main(x, 'gain'),
         main(x, 'loss')
      )
   }
}

#calcChromArmCnChange(df[1,grep('^ploidy',colnames(df))])

