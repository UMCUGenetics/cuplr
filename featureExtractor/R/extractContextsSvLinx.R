#' Extract SV contexts from LINX
#'
#' @param vis.sv.data Path to the *vis_sv_data.tsv file or a dataframe of this file
#'
#' @return An integer vector of counts
#' @export
#'
extractContextsSvLinx <- function(vis.sv.data=NULL){
   
   ## Load inputs --------------------------------
   if(is.character(vis.sv.data)){
      linx_svs <- read.delim(vis.sv.data, check.names=F)
   } else if(is.data.frame(vis.sv.data)){
      linx_svs <- vis.sv.data
      rm(vis.sv.data)
   } else {
      stop('`vis.sv.data` must be a path or a dataframe')
   }
   
   ## Main --------------------------------
   ## DEL, DUP and COMPLEX SVs stratified by length
   DEL_DUP <- countDelDupByLen(linx_svs)[1,]
   COMPLEX <- countComplexByLen(linx_svs)[1,]
   
   ## Number of LINEs
   LINEs <- length(unique(
      subset(linx_svs, ResolvedType=='LINE', ClusterId, drop=T)
   ))
   
   ## Number of double minutes
   double_minutes <- length(unique(
      subset(linx_svs, ResolvedType=='DOUBLE_MINUTE', ClusterId, drop=T)
   ))
   
   ## Total SV load
   sv_load <- length(unique(linx_svs$ClusterId))
   
   ## 
   COMPLEX.largest_cluster <- (function(){
      df <- subset(linx_svs, ResolvedType=='COMPLEX')
      if(nrow(df)==0){ return(0) }
      cluster_sizes <- table(df$ClusterId)
      max(cluster_sizes)
   })()
   
   ## Output --------------------------------
   out <- c(
      sv_load=sv_load,
      LINEs=LINEs,
      DEL_DUP,
      COMPLEX.largest_cluster=COMPLEX.largest_cluster,
      COMPLEX,
      double_minutes=double_minutes
   )
   out[is.na(out)] <- 0
   return(out)
}

####################################################################################################
#' Count DELs and DUPs stratified by length
#'
#' @param linx.svs A dataframe of a *linx.vis_sv_data.tsv file
#' @param bin.breaks An integer vector specifying the length bin breaks
#' @param sample.id.colname A string specifying the column name used by indicating the sample name. 
#' If unspecified, all SVs are assumed to be from the same sample.
#' @param output 'contexts' returns a matrix of the counts of each context. 'raw' returns a 
#' dataframe with annotations for each SV (i.e. row).
#'
#' @return An integer vector
#' @export
#'
countDelDupByLen <- function(
   linx.svs, bin.breaks=c(0,300,10^(3:7),Inf), 
   sample.id.colname=NULL, output='contexts'
){
   
   if(!(output %in% c('contexts','raw'))){
      stop("`output` must be 'contexts' or 'raw'")
   }
   
   ## Template context output --------------------------------
   bin_names <- levels( cut(0L, bin.breaks, right=FALSE, include.lowest=FALSE) )
   context_names <-  c(
      paste0('DEL_', bin_names),
      paste0('DUP_', bin_names)
   )
   context_counts_template <- structure(
      rep(0, length(context_names)), 
      names=context_names
   )
   context_counts_template <- t(context_counts_template)
   
   ## Prep data --------------------------------
   df <- linx.svs

   ##
   if(nrow(df)==0){ return(context_counts_template) }
   if(nrow(df)==0 & output=='raw'){ stop('`linx.svs` cannot have no rows when output="raw"') }
   
   ## Get or make sample column
   if(is.null(sample.id.colname)){
      df$sample <- 'sample1'
   } else {
      df$sample <- df[[sample.id.colname]]
      #df[[sample.id.colname]] <- NULL
   }
   df$sample <- as.factor(df$sample)
   
   ## Subset for relevant rows and cols
   df <- df[
      df$ResolvedType %in% c('DEL','DUP'),
      c('sample', 'ClusterId', 'ResolvedType', 'ChrStart', 'ChrEnd', 'PosStart', 'PosEnd')
   ]
   
   if(nrow(df)==0){ return(context_counts_template) }
   
   ## Use sample specific cluster ids
   df$ClusterId <- paste0( as.integer(df$sample),'_',df$ClusterId )
   
   ## Main --------------------------------
   ## Select intrachromosomal SVs
   df$is_same_chrom <- df$ChrStart==df$ChrEnd
   
   ## All events in cluster are on the same chromosome?
   is_intrachrom_cluster <- sapply(split(df$is_same_chrom, df$ClusterId), all)
   df$is_intrachrom_cluster <- is_intrachrom_cluster[ df$ClusterId ]
   df <- df[df$is_intrachrom_cluster,]
   if(nrow(df)==0){ return(context_counts_template) }
   
   # intrachrom_clusters <- aggregate(
   #    df$is_same_chrom, 
   #    list(ClusterId=df$ClusterId), 
   #    all
   # ) ## All events in cluster are on the same chromosome?
   # 
   # df <- df[
   #    df$ClusterId %in% intrachrom_clusters$ClusterId[ is_intrachrom_cluster ]
   # ,]
   
   NULL -> df$ChrStart -> df$ChrEnd -> df$is_same_chrom -> df$is_intrachrom_cluster
   
   ## Flatten clusters of SVs into one row
   df$is_start_fragment <- !duplicated(df$ClusterId, fromLast=F)
   df$is_end_fragment <- !duplicated(df$ClusterId, fromLast=T)
   
   df_flat <- data.frame(
      sample=df$sample[df$is_start_fragment],
      ClusterId=df$ClusterId[df$is_start_fragment],
      ResolvedType=df$ResolvedType[df$is_start_fragment],
      start=df$PosStart[df$is_start_fragment],
      end=df$PosEnd[df$is_end_fragment]
   )
   
   ## Calculate SV lengths
   # ## start coords are 1-based
   # all(df_flat$start <= df_flat$end)
   # df_flat[df_flat$start == df_flat$end,]
   df_flat$length <- 1 + df_flat$end - df_flat$start
   
   ## Bin SVs by length
   df_flat$length_bin <- cut(df_flat$length, bin.breaks, right=FALSE, include.lowest=FALSE)
   
   ## Remove SVs not falling into the bins specified in bin.breaks
   df_flat <- df_flat[!is.na(df_flat$length_bin),]
   if(nrow(df_flat)==0){ return(context_counts_template) }
   
   ## Make SV type/length contexts
   df_flat$context <- paste0(df_flat$ResolvedType, '_', df_flat$length_bin)
   df_flat <- df_flat[order(df_flat$ResolvedType, df_flat$length_bin),]
   df_flat$context <- factor(df_flat$context, context_names)
   
   if(output=='raw'){ return(df_flat) }
   
   ## Make context matrix -------------------------------
   m <- table(df_flat$sample, df_flat$context)
   m <- unclass(m)
   
   return(m)
}

####################################################################################################
#' Count COMPLEX SVs stratified by length
#'
#' @param linx.svs A dataframe of a *linx.vis_sv_data.tsv file
#' @param bin.breaks An integer vector specifying the length bin breaks.
#' @param sample.id.colname A string specifying the column name used by indicating the sample name. 
#' If unspecified, all SVs are assumed to be from the same sample.
#' @param output 'contexts' returns a matrix of the counts of each context. 'raw' returns a 
#' dataframe with annotations for each cluster.
#'
#' @return An integer vector
#' @export
#'
countComplexByLen <- function(
   linx.svs, bin.breaks=c(0,25,50,100,200,400,800,Inf), 
   sample.id.colname=NULL, output='contexts'
){
   
   if(!(output %in% c('contexts','raw'))){
      stop("`output` must be 'contexts' or 'raw'")
   }
   
   ## Template context output --------------------------------
   bin_names <- levels( cut(0L, bin.breaks, right=FALSE, include.lowest=FALSE) )
   context_names <-  paste0('COMPLEX_', bin_names)
   context_counts_template <- structure(
      rep(0, length(context_names)), 
      names=context_names
   )
   context_counts_template <- t(context_counts_template)
   
   ## Prep data --------------------------------
   df <- linx.svs
   
   ##
   if(nrow(df)==0){ return(context_counts_template) }
   if(nrow(df)==0 & output=='raw'){ stop('`linx.svs` cannot have no rows when output="raw"') }
   
   ## Get or make sample column
   if(is.null(sample.id.colname)){
      df$sample <- 'sample1'
   } else {
      df$sample <- df[[sample.id.colname]]
   }
   df$sample <- as.factor(df$sample)
   
   ## Subset for COMPLEX SVs
   df <- subset(df, ResolvedType=='COMPLEX', c(sample, ClusterId, SvId))
   
   if(nrow(df)==0){ return(context_counts_template) }
   
   ## Bin complex clusters by length --------------------------------
   complex_len <- with(df,{
      out <- aggregate(
         ClusterId,
         list(sample=sample, ClusterId=ClusterId),
         function(x){ length(x) }
      )
      colnames(out)[length(out)] <- 'len'
      return(out)
   })
   
   complex_len$len_bin <- cut(complex_len$len, bin.breaks, right=FALSE, include.lowest=FALSE)

   if(output=='raw'){ return(complex_len) }
   
   ## Make context matrix -------------------------------
   m <- table(complex_len$sample, complex_len$len_bin)
   m <- unclass(m)
   colnames(m) <- paste0('COMPLEX_',colnames(m))
   return(m)
}







