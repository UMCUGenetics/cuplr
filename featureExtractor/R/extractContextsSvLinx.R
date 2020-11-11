#' Make interval lookup table
#'
#' @param bins A numeric vector specifying the interval breaks
#'
#' @return A dataframe
#' @export
#'
mkIntervalLookup <- function(bins=c(0, 10^c(3:7))){
   intervals <- data.frame(
      lower=bins[-length(bins)],
      upper=bins[-1]
   )
   
   intervals$name <- paste0(
      formatC(intervals$lower, format='e', digits=0),
      '_',
      formatC(intervals$upper, format='e', digits=0),
      '_bp'
   )
   
   intervals$name <- gsub('[+]','',intervals$name)
   intervals$counts <- 0
   return(intervals)
}

#' Count binned SV length
#'
#' @param v.sv.len A numeric vector of SV lengths
#' @param bins A numeric vector specifying the interval breaks
#'
#' @return An integer vector of counts
#' @export
#'
countBinnedSvLen <- function(v.sv.len, bins=c(0, 10^c(3:7), Inf)){
   #v.sv.len=df[df$ResolvedType=='DEL','sv_len']
   #bins=len.bins.del
   
   intervals <- mkIntervalLookup(bins)
   
   if(length(v.sv.len)!=0){
      intervals$counts <- unlist(lapply(1:nrow(intervals), function(i){
         sum(intervals$lower[i] <= v.sv.len & v.sv.len < intervals$upper[i])
      }))
   }
   
   structure(
      intervals$counts,
      names=intervals$name
   )
}

####################################################################################################
#' Extract SV contexts from LINX
#'
#' @param vis.sv.data.path Path to the *vis_sv_data.tsv file
#' @param df A dataframe of the *vis_sv_data.tsv file
#' @param len.bins.del A numeric vector specifying the interval breaks to bin DEL lengths
#' @param len.bins.dup A numeric vector specifying the interval breaks to bin DUP lengths
#' @param resolved.type.annotations.path Path to the txt file containing annotations for each LINX
#' ResolvedType
#' @param sel.cols A character vector with the names: ResolvedType, ClusterId, PosStart, PosEnd.
#' The value corresponding to each name should refer to a column name in the txt file. This is used 
#' to translate the column names in the txt file to the column names that the function will use.
#'
#' @return An integer vector of counts
#' @export
#'
extractContextsSvLinx <- function(
   vis.sv.data.path=NULL, df=NULL, 
   len.bins.del=c(0, 10^c(3:7), Inf), len.bins.dup=len.bins.del,
   resolved.type.annotations.path=RESOLVED_TYPE_ANNOTATIONS,
   sel.cols=NULL
){
   #vis.sv.data.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/181114_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.vis_sv_data.tsv'
   #vis.sv.data.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/160805_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.vis_sv_data.tsv'
   #vis.sv.data.path="/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/featureExtractor/test/COLO829/raw//COLO829v003T.linx.vis_sv_data.tsv"
   #df=vis_sv_data[vis_sv_data$SampleId=='XXXXXXXX',]
   #len.bins.del=c(0, 10^c(3:7), Inf)
   #len.bins.dup=len.bins.del
   #resolved.type.annotations.path=RESOLVED_TYPE_ANNOTATIONS
   
   ## Load inputs --------------------------------
   if(!is.null(vis.sv.data.path)){ df <- read.delim(vis.sv.data.path, stringsAsFactors=F) }
   
   required_cols <- c('ResolvedType','ClusterId','PosStart','PosEnd')
   df <- selectRequiredCols(df, required.cols=required_cols, sel.cols=sel.cols)

   resolved_type_annotations <- read.delim(resolved.type.annotations.path, stringsAsFactors=F)
   
   ## Initialize output --------------------------------
   interval_names_del <- paste0('DEL_',mkIntervalLookup(len.bins.del)$name)
   interval_names_dup <- paste0('DUP_',mkIntervalLookup(len.bins.dup)$name)
   
   other_feature_names <- resolved_type_annotations$ResolvedType
   other_feature_names <- other_feature_names[!(other_feature_names %in% c('DEL','DUP','COMPLEX'))]
      
   feature_names <- c(
      interval_names_del,
      interval_names_dup,
      'COMPLEX',
      other_feature_names
   )
   
   counts <- structure(
      rep(0, length(feature_names)),
      names=feature_names
   )
   
   ## DEL/DUP length --------------------------------
   df$sv_len <- abs(df$PosEnd - df$PosStart)
   sv_type_len <- list()
   
   sv_type_len$DEL <- countBinnedSvLen(
      df[df$ResolvedType=='DEL','sv_len'],
      bins=len.bins.del
   )

   sv_type_len$DUP <- countBinnedSvLen(
      df[df$ResolvedType=='DUP','sv_len'],
      bins=len.bins.del
   )
   
   sv_type_len <- do.call(c, sv_type_len)
   names(sv_type_len) <- gsub('[.]','_',names(sv_type_len))
   
   counts[names(sv_type_len)] <- sv_type_len
   
   ## Count complex events --------------------------------
   df_complex <- df[df$ResolvedType=='COMPLEX',]
   counts['COMPLEX'] <- length(unique(df_complex[,'ClusterId'])) ## Consider each event cluster as one event
   

   ## Other resolved types --------------------------------
   df_other <- df[
      df$ResolvedType %in% other_feature_names,
      c('ClusterId','ResolvedType')
   ]
   
   df_other <- unique(df_other)
   tab <- table(df_other$ResolvedType)
   
   counts[names(tab)] <- as.numeric(tab)
   
   ## Return --------------------------------
   counts <- counts[feature_names]
   return(counts)
}








