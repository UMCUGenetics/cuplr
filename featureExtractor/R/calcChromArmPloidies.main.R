#' Calculate overall chrom arm copy numbers
#'
#' @description This function first rounds copy numbers (CN) to integers so that CN segments can be
#' grouped together. Per chrom arm, the coverage of each CN category is calculated (i.e.
#' cumulative segment size). The chrom arm CN is (roughly) defined as the CN category with the
#' highest cumulative segment size
#'
#' @param cnv.file Path to purple cnv file
#' @param cnv purple cnv file loaded as a dataframe
#' @param out.file Path to output file. If NULL, returns a named vector
#' @param sel.cols A character vector with the names: chrom, start, end, total_cn, major_cn,
#' minor_cn. The value corresponding to each name should refer to a column name in the txt file.
#' This is used to translate the column names in the txt file to the column names that the function
#' will use.
#' @param mode Can be total_cn (for determining arm CN changes) or minor_cn (for determining arm LOH)
#' @param min.rel.cum.segment.size If a chrom arm has a CN category that covers >0.5 (i.e 50%;
#' default) of a chrom arm, this CN is the copy number of the arm
#' @param max.rel.cum.segment.size.diff This value (default 0.1) determines whether which CN
#' categories are considered to cover equal lengths of the chrom arm. For example, (by default) 2
#' CN categories covering 0.40 and 0.31 of a chrom arm are considered equally contributing. When
#' these CN categories have similar cumulative segment size as the one with the highest, if one of
#' these have the same CN as the genome CN, return the genome CN. Otherwise, simply return the one
#' with the highest segment support (as is done above).
#' @param keep.chroms Which chromosomes to keep?
#' @param one.armed.chroms Which chromosomes are considered to have only one arm?
#' @param verbose Show progress messages?
#'
#' @return A named vector of chrom arm copy numbers, or specified writes a table to out.file if
#' specified
#' @export
#'
calcChromArmPloidies <- function(
   ## I/O
   cnv.file=NULL, cnv=NULL, out.file=NULL, sel.cols=NULL,

   ## Params
   mode=c('total_cn','minor_cn'),
   min.rel.cum.segment.size=if(mode[1L]=='minor_cn'){ 0.9 } else { 0.5 },
   max.rel.cum.segment.size.diff=0.1,

   ## Misc
   keep.chroms=paste0('chr',c(1:22,'X')),
   one.armed.chroms=c(13,14,15,21,22),
   verbose=F
){

   ## Load data --------------------------------
   if(!is.null(cnv.file)){
      cnv <- read.delim(cnv.file, check.names=F, stringsAsFactors=F)
   }

   if(is.null(sel.cols)){
      sel.cols <- c(
         chrom='chromosome',start='start',end='end',
         total_cn='copyNumber',major_cn='majorAllelePloidy',minor_cn='minorAllelePloidy'
      )
   }

   cnv <- selectRequiredCols(
      df=cnv,
      required.cols=c('chrom','start','end','total_cn','major_cn','minor_cn'),
      sel.cols=sel.cols
   )

   GenomeInfoDb::seqlevelsStyle(cnv$chrom)<- 'NCBI'
   if(!is.null(keep.chroms)){
      GenomeInfoDb::seqlevelsStyle(keep.chroms)<- 'NCBI'
      cnv <- cnv[cnv$chrom %in% keep.chroms,]
   }

   ## Pre-calculations --------------------------------
   cnv$segment_size <- (cnv$end - cnv$start) + 1

   cnv$total_cn[cnv$total_cn < 0] <- 0 ## Ceiling negative values
   cnv$total_cn_int <- round(cnv$total_cn)

   cnv$minor_cn[cnv$minor_cn < 0] <- 0
   cnv$minor_cn_int <- round(cnv$minor_cn)

   ##----------------------------------------------------------------
   if(verbose){ message('Splitting CN segments by chrom arm...') }
   cnv$chrom_arm <- getChromArm(
      data.frame(chrom=cnv$chrom, start=cnv$start, end=cnv$end),
      one.armed.chroms=one.armed.chroms,
      centro.intervals.rough.fix=T
   )
   cnv$chrom_arm <- gsub('chr','',cnv$chrom_arm)
   cnv$chrom_arm <- factor(cnv$chrom_arm, unique(cnv$chrom_arm))
   cnv_split <- split(cnv, cnv$chrom_arm)

   minor_cn_segment_support <- lapply(cnv_split, function(i){
      #i=cnv_split$`Yp`
      df <- aggregate(i$segment_size, by=list(i$minor_cn_int), FUN=sum)
      colnames(df) <- c('minor_cn_int','cum_segment_size')
      df <- as.data.frame(lapply(df, as.integer)) ## aggregate can return lists instead of vectors as output

      df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)

      #df[which.max(df$cum_segment_size_rel),]
      if(nrow(df)!=0L){
         df <- df[order(df$cum_segment_size_rel, decreasing=T),]
         return(df)
      } else {
         df[1L,] <- NA
         return(df)
      }
   })

   ##----------------------------------------------------------------
   if(verbose){ message('Calculating copy number segment support...') }
   calcSegmentSupport <- function(colname){
      #colname='total_cn_int'
      lapply(cnv_split, function(i){
         #i=cnv_split$`1q`
         df <- aggregate(i$segment_size, by=list(i[,colname]), FUN=sum)
         colnames(df) <- c('cn','cum_segment_size')
         #df <- as.data.frame(lapply(df, as.integer)) ## aggregate can return lists instead of vectors as output

         df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)

         #df[which.max(df$cum_segment_size_rel),]
         if(nrow(df)!=0L){
            df <- df[order(df$cum_segment_size_rel, decreasing=T),]
            return(df)
         } else {
            df[1L,] <- NA
            return(df)
         }
      })
   }

   if(mode[1L]=='total_cn'){
      cn_segment_support <- calcSegmentSupport('total_cn_int')
   } else if(mode[1L]=='minor_cn'){
      cn_segment_support <- calcSegmentSupport('minor_cn_int')
   } else {
      stop("`mode` must be 'total_cn' or 'minor_cn'")
   }

   ## CN with most frequent total segment support is preliminary CN
   arm_cn_prelim <- unlist(lapply(cn_segment_support, function(i){ i[1L,'cn'] }))

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

         if(is.na(i[1L,1L])){ return(NA) }

         ## E.g. >=50% of arm has CN of 2, then this is the CN
         if( i[1L,'cum_segment_size_rel'] >= min.rel.cum.segment.size ){
            return(i[1L,'cn'])
         }

         ## When multiple CNs have similar segment support as the one with the highest, if one
         ## of these have the same CN as the genome CN, return the genome CN. Otherwise, simply
         ## return the one with the highest segment support (as is done above)
         i$diffs <- i[1L,'cum_segment_size_rel'] - i[,'cum_segment_size_rel']
         cn_doubt <- i[i$diffs < max.rel.cum.segment.size.diff,'cn']

         if(any(cn_doubt==genome_cn)){
            return(genome_cn)
         } else {
            return(i[1L,'cn'])
         }
      }))
   }

   ploidies <- c(arm_cn, genome=genome_cn)

   ## Ensure consistent output (e.g. for when CN data is missing) --------------------------------
   if(verbose){ message('Preparing output...') }
   chrom_arm_names <- paste0(rep(keep.chroms,each=2L),c('p','q'))
   chrom_arm_names <- chrom_arm_names[!(chrom_arm_names %in% paste0(one.armed.chroms,'p'))] ## Keep only q arm for one arm chromosomes
   chrom_arm_names <- c(chrom_arm_names,'genome')

   ## Output --------------------------------
   out <- structure(rep(NA,length(chrom_arm_names)), names=chrom_arm_names)
   out[names(ploidies)] <- ploidies
   
   if(is.null(out.file)){
      return(out)
   } else {
      out <- data.frame(chrom=names(out),ploidy=out,row.names=NULL)
      if(verbose){ message('Writing output...') }
      write.table(
         out,
         if(grepl('[.]gz$',out.file)){ gzfile(out.file) } else { out.file },
         sep='\t', quote=F, row.names=F,
      )
   }

}

################################################################################
#' Calculate change in chrom arm copy number compared to genome copy number
#'
#' @rdname calcChromArmCnChange
#'
#' @param x Output from calcChromArmPloidies() as a matrix or vector
#' @param direction Can be'gain' or 'loss'
#'
#' @return A matrix or vector
#' @export
#'
calcChromArmCnChange <- function (x, ...) {
   UseMethod("calcChromArmCnChange", x)
}

#' @rdname calcChromArmCnChange
#' @method calcChromArmCnChange numeric
#' @export
calcChromArmCnChange.numeric <- function(x, direction){
   #x <- unlist(features$ploidy[1,])

   genome_index <- which(names(x)=='genome')
   arm <- x[-genome_index]
   genome <- x[genome_index]

   x_diff <- arm - genome
   if(direction=='gain'){
      x_diff[x_diff < 0] <- 0
   } else if(direction=='loss') {
      x_diff[x_diff > 0] <- 0
   } else {
      stop("`direction` must be 'gain' or 'loss'")
   }
   x_diff <- abs(x_diff)

   return(x_diff)
}

#' @rdname calcChromArmCnChange
#' @method calcChromArmCnChange matrix
#' @export
calcChromArmCnChange.matrix <- function(m, direction){
   #m=features$ploidy

   genome_index <- which(colnames(m)=='genome')
   arm <- m[,-ncol(m)]
   genome <- m[,ncol(m)]

   x_diff <- arm - genome
   if(direction=='gain'){
      x_diff[x_diff < 0] <- 0
   } else if(direction=='loss') {
      x_diff[x_diff > 0] <- 0
   } else {
      stop("`direction` must be 'gain' or 'loss'")
   }
   x_diff <- abs(x_diff)

   return(x_diff)
}

#' @rdname calcChromArmCnChange
#' @method calcChromArmCnChange data.frame
#' @export
calcChromArmCnChange.data.frame <- calcChromArmCnChange.matrix



