#' Assign clusters to a vector of values based on run length encoding
#'
#' @param x An integer vector
#'
#' @return
#' @export
#'
rle2Clusters <- function(x){
   rle_out <- rle(x)
   unlist(lapply(1:length(rle_out$lengths), function(i){
      rep(i, rle_out$lengths[i])
   }))
}

####################################################################################################
#' Converts chrom/pos into linear coordinates
#'
#' @param chrom A character vector indicating the chromosome
#' @param pos An integer vector indicating the position
#' @param ref.genome A BSgenome object
#'
#' @return An integer vector of linear coordinates
#' @export
#'
linearizeChromPos <- function(chrom, pos, ref.genome=mutSigExtractor::DEFAULT_GENOME){
   #chrom=kat_snv$chrom
   #pos=kat_snv$pos
   
   if(length(chrom)!=length(pos)){ stop('`chrom` and `pos` are not the same length') }
   
   chrom <- as.character(chrom)
   GenomeInfoDb::seqlevelsStyle(chrom)<-  GenomeInfoDb::seqlevelsStyle(chrom)
   
   chrom_lengths <- (function(){
      v <- GenomeInfoDb::seqlengths(ref.genome)
      v[names(v) %in% chrom]
   })()
   
   offsets <- cumsum(as.numeric(chrom_lengths))
   offsets <- c(0,offsets[-1])
   names(offsets) <- names(chrom_lengths)
   
   unlist(Map(function(chrom,pos){
      pos + offsets[chrom]
   }, chrom, pos), use.names=F)
}

####################################################################################################
#' Combine feature columns by summing them together
#'
#' @param x A matrix or dataframe 
#' @param target.features A character vector of column (feature) names
#' @param regex Instead of target.features, a regex can be specified
#' @param target.name A name to assign the new feature. If unspecified, will use the name of the 
#' first old feature
#'
#' @return The original matrix or dataframe with the indicated features combined
#' @export
#'
combineFeatures <- function(x, target.features=NULL, regex=NULL, target.name=NULL){
   
   if(!is.null(regex)){
      target.features <- grep(regex,colnames(x),value=T)
      if(length(target.features)==0){ stop('No features match the regex pattern') }
   }
   
   target_col <- target.features[1]
   
   x[,target_col] <- rowSums(x[,target.features])
   
   rm_cols <- target.features[target.features!=target_col]
   x <- x[,!(colnames(x) %in% rm_cols)]
   
   if(!is.null(target.name)){
      colnames(x)[colnames(x)==target_col] <- target.name
   }
   
   return(m)
}