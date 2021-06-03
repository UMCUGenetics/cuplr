#' Extract SV contexts from LINX
#'
#' @param vis.sv.data Path to the *vis_sv_data.tsv file or a dataframe of this file
#'
#' @return An integer vector of counts
#' @export
#'
extractContextsSvLinx <- function(vis.sv.data=NULL){
   # if(F){
   #    vis.sv.data="/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/DO1003-from-jar/linx14/DO1003T.linx.vis_sv_data.tsv"
   #    vis.sv.data="/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104-update3//somatics/161205_HMFregXXXXXXXX/linx14/XXXXXXXX.linx.vis_sv_data.tsv"
   #    vis.sv.data="/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/160601_HMFreg0056_FR10302053_FR10302054_XXXXXXXX/linx14/XXXXXXXX.linx.vis_sv_data.tsv"
   # }
   
   ## Load inputs --------------------------------
   if(is.character(vis.sv.data)){
      linx_svs <- read.delim(vis.sv.data, check.names=F)
   } else if(is.data.frame(vis.sv.data)){
      linx_svs <- vis.sv.data
      rm(vis.sv.data)
   } else {
      stop('`vis.sv.data` must be a path or a dataframe')
   }
   
   #linx_svs$IsSynthetic <- as.logical(linx_svs$IsSynthetic)
   linx_svs$InDoubleMinute <- as.logical(linx_svs$InDoubleMinute)
   linx_svs$HasFoldback <- with(linx_svs, InfoStart=='FOLDBACK' | InfoEnd=='FOLDBACK' | ResolvedType=='RESOLVED_FOLDBACK' )
   
   ## Punctuated structural events (use raw counts) --------------------------------
   foldbacks <- length(unique(
      subset(linx_svs, HasFoldback, ClusterId, drop=T)
   ))
   
   double_minutes <- length(unique(
      subset(linx_svs, InDoubleMinute, ClusterId, drop=T)
   ))
   
   svs_complex <- subset(linx_svs, ResolvedType=='COMPLEX')
   if(nrow(svs_complex)!=0){
      complex_lengths <- with(subset(linx_svs, ResolvedType=='COMPLEX'),{
         agg <- aggregate(
            ClusterId,
            list(ClusterId=ClusterId),
            function(x){ length(x) }
         )
         colnames(agg)[ncol(agg)] <- 'length'
         return(agg)
      })
      
      complex.n_events <- length(complex_lengths$ClusterId)
      complex.largest_cluster <- max(complex_lengths$length)
   } else {
      complex.n_events <- 0
      complex.largest_cluster <- 0
   }
   
   ## Ongoing structural events (use relative counts) --------------------------------
   sv_load <- length(unique(linx_svs$ClusterId))
   
   LINEs <- length(unique(
      subset(linx_svs, ResolvedType=='LINE', ClusterId, drop=T)
   ))
   
   simple_events <- countSimpleSvEvents(linx_svs)
   
   ## Output --------------------------------
   out <- c(
      n_events=sv_load,
      complex.n_events=complex.n_events,
      complex.largest_cluster=complex.largest_cluster,
      double_minutes=double_minutes,
      foldbacks=foldbacks,
      LINEs=LINEs/sv_load,
      simple_events/sv_load
   )
   out[is.na(out)] <- 0
   return(out)
}

####################################################################################################
#' Count DELs and DUPs stratified by length
#'
#' @param linx.svs A dataframe of a *linx.vis_sv_data.tsv file
#' @param bin.breaks An integer vector specifying the length bin breaks
#'
#' @return An integer vector
#' @export
#'
countSimpleSvEvents <- function(linx.svs, bin.breaks=c(0,10^(3:7),Inf)){
   
   bin_intervals <- levels(cut(0, bin.breaks, right=FALSE, include.lowest=FALSE))
   context_names <- c(
      paste0('DEL_',bin_intervals),
      paste0('DUP_',bin_intervals)
   )
   
   context_counts <- structure(
      rep(0,length(context_names)),
      names=context_names
   )
   
   df <- subset(
      linx.svs, 
      ResolvedType %in% c('DEL','DUP'), 
      c(ClusterId, ResolvedType, ChrStart, ChrEnd, PosStart, PosEnd)
   )
   
   if(nrow(df)==0){ return(context_counts) }
   
   ## Select intrachromosomal SVs
   df$is_same_chrom <- df$ChrStart==df$ChrEnd
   cluster_in_same_chrom <- aggregate(df$is_same_chrom, list(ClusterId=df$ClusterId), all) ## All events in cluster are on the same chromosome?
   
   df <- df[
      df$ClusterId %in% cluster_in_same_chrom$ClusterId[ cluster_in_same_chrom$x ]
      ,]
   
   NULL -> df$ChrStart -> df$ChrEnd -> df$is_same_chrom -> cluster_in_same_chrom
   
   ## Calculate SV lengths
   pos_min <- aggregate(df$PosStart, list(ClusterId=df$ClusterId), min)
   pos_max <- aggregate(df$PosEnd, list(ClusterId=df$ClusterId), max)
   
   df_flat <- unique( df[,c('ClusterId', 'ResolvedType')] )
   df_flat$start <- pos_min$x[ match(df_flat$ClusterId, pos_min$ClusterId) ]
   df_flat$end <- pos_max$x[ match(df_flat$ClusterId, pos_max$ClusterId) ]
   df_flat$length <- 1 + df_flat$end - df_flat$start
   ## all(df_flat$start <= df_flat$end)
   ## df_flat[df_flat$start == df_flat$end,]
   ## start coords are 1-based
   
   ## Make SV type/length contexts
   df_flat$length_bin <- cut(df_flat$length, bin.breaks, right=FALSE, include.lowest=FALSE)
   df_flat$context <- paste0(df_flat$ResolvedType, '_', df_flat$length_bin)
   
   ## Make context vector
   tab <- table(df_flat$context)
   context_counts[names(tab)] <- tab
   return(context_counts)
}







