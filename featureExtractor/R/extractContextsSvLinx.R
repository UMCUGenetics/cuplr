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
#' Assign COMPLEX SV subtypes
#'
#' @param df A dataframe of the *vis_sv_data.tsv file
#' @param verbose Show progress?
#'
#' @return A data.frame or character vector
#' @export
#'
annotateComplexSVs <- function(df, verbose=F){
   #df=df_complex[df_complex$SampleId=='XXXXXXXX',]
   
   df$Ploidy[df$Ploidy < 0 ] <- 0
   df$Ploidy <- round(df$Ploidy)
   
   df_split <- split(df, factor(df$ClusterId, unique(df$ClusterId)))
   if(verbose){ 
      message('Calculating stats for each cluster...') 
      pb <- txtProgressBar(max=length(df_split), style=3)
      counter <- 0
   }
   
   cluster_ann <- do.call(rbind, lapply(df_split, function(i){
      #i=df_split[[2]]
      
      if(verbose){ 
         counter <<- counter + 1
         setTxtProgressBar(pb,counter) 
      }
      
      n_sv <- nrow(i)
      
      ## Chrom
      tab_chrom <- sort(table(i$ChrStart), decreasing=T)
      n_chrom <- length(tab_chrom)
      chrom_prop <- structure(
         round(tab_chrom[1:3] / n_sv, 3),
         names=paste0('chrom_',1:3,'_prop')
      )
      
      ## CN states
      tab_cn <- sort(table(i$Ploidy), decreasing=T)
      n_cn_states <- length(tab_cn)
      cn_state_prop <- structure(
         round(tab_cn[1:3] / n_sv, 3),
         names=paste0('cn_state_',1:3,'_prop')
      )
      
      foldback_prop <- round(sum(i$InfoStart=='FOLDBACK') / n_sv, 3)
      max_cn <- max(i$Ploidy)
      
      c(
         n_sv=n_sv,
         foldback_prop=foldback_prop, 
         
         n_chrom=n_chrom, 
         chrom_prop,
         
         max_cn=max_cn, 
         n_cn_states=n_cn_states,
         cn_state_prop
      )
   }))
   cluster_ann[is.na(cluster_ann)] <- 0
   
   cluster_ann <- cbind(
      ClusterId=rownames(cluster_ann),
      as.data.frame(cluster_ann)
   )
   rownames(cluster_ann) <- NULL
   
   if(verbose){ message('\nDetermining complex type...') }
   complex_type <- data.frame(
      foldbacks = with(cluster_ann,{ 
         foldback_prop>=0.6 &
            n_sv>=5
      }),
      
      hiSvLoad = cluster_ann$n_sv>=50,
      
      fewChrom_ploidyOscill = with(cluster_ann,{ 
         n_sv>=10 &
            n_chrom<=2 & 
            n_cn_states>=2 & 
            (cn_state_1_prop + cn_state_2_prop)>=0.8
      }),
      
      multiChrom_ploidyStable = with(cluster_ann,{ 
         n_chrom>=3 & 
            cn_state_1_prop>=0.8 & 
            chrom_1_prop>=0.2 & chrom_2_prop>=0.2
      }),
      
      ploidyUnstable = cluster_ann$n_cn_states >= 10
   )
   complex_type$other <- !apply(complex_type, 1, any)
   
   cluster_ann$type <- colnames(complex_type)[ max.col(complex_type) ]
   cluster_ann$type <- factor(cluster_ann$type, colnames(complex_type))
   
   return(cluster_ann)
}

COMPLEX_SVTYPES <- c('foldbacks','hiSvLoad','fewChrom_ploidyOscill','multiChrom_ploidyStable','ploidyUnstable','other')

####################################################################################################
#' Extract SV contexts from LINX
#'
#' @param vis.sv.data.path Path to the *vis_sv_data.tsv file
#' @param df A dataframe of the *vis_sv_data.tsv file
#' @param len.bins.del A numeric vector specifying the interval breaks to bin DEL lengths
#' @param len.bins.dup A numeric vector specifying the interval breaks to bin DUP lengths
#' @param resolved.type.annotations.path Path to the txt file containing annotations for each LINX
#' ResolvedType
#' @param split.complex If TRUE, COMPLEX SVs will be split into foldbacks, hiSvLoad, 
#' fewChrom_ploidyOscill, multiChrom_ploidyStable, ploidyUnstable, and other
#'
#' @return An integer vector of counts
#' @export
#'
extractContextsSvLinx <- function(
   vis.sv.data.path=NULL, df=NULL, 
   len.bins.del=c(0, 10^c(3:7), Inf), len.bins.dup=len.bins.del,
   resolved.type.annotations.path=RESOLVED_TYPE_ANNOTATIONS,
   split.complex=F
){
   #vis.sv.data.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/181114_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.vis_sv_data.tsv'
   #vis.sv.data.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/160805_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.vis_sv_data.tsv'
   #df=vis_sv_data[vis_sv_data$SampleId=='XXXXXXXX',]
   #len.bins.del=c(0, 10^c(3:7), Inf)
   #len.bins.dup=len.bins.del
   #resolved.type.annotations.path=RESOLVED_TYPE_ANNOTATIONS
   
   if(!is.null(vis.sv.data.path)){ df <- read.delim(vis.sv.data.path, stringsAsFactors=F) }
   
   resolved_type_annotations <- read.delim(resolved.type.annotations.path, stringsAsFactors=F)
   
   ## Initialize output --------------------------------
   interval_names_del <- paste0('DEL_',mkIntervalLookup(len.bins.del)$name)
   interval_names_dup <- paste0('DUP_',mkIntervalLookup(len.bins.dup)$name)
   
   other_feature_names <- resolved_type_annotations$ResolvedType
   other_feature_names <- other_feature_names[!(other_feature_names %in% c('DEL','DUP','COMPLEX'))]
      
   feature_names <- c(
      interval_names_del,
      interval_names_dup,
      if(split.complex){ paste0('COMPLEX_',COMPLEX_SVTYPES) } else { 'COMPLEX' },
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
   
   ## Split complex --------------------------------
   complex_types <- table( annotateComplexSVs(df)$type )
   names(complex_types) <- paste0('COMPLEX_',names(complex_types))

   if(split.complex){
      counts[names(complex_types)] <- complex_types
   } else {
      counts['COMPLEX'] <- sum(complex_types)
   }

   ## Other resolved types --------------------------------
   df_ss <- df[
      df$ResolvedType %in% other_feature_names,
      c('ClusterId','ResolvedType')
   ]
   
   df_ss <- unique(df_ss)
   tab <- table(df_ss$ResolvedType)
   
   counts[names(tab)] <- as.numeric(tab)
   
   ## Return --------------------------------
   counts <- counts[feature_names]
   return(counts)
}








