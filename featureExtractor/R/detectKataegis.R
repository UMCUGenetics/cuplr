#' Detect kataegis from SNVs
#' 
#' @description Detects kataegis from SNVs provided from a vcf file or dataframe. A piecewise 
#' constant fit (using rpart (R in built regression tree algorithm)) is performed on the
#' intermutation distances. By default, segments with min.muts>=6 and imd.thres<=1000 are considered
#' kataegis foci. These default values were chosen according to d'Antonio et al Cell Rep. 2016 
#' (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4972030/)
#' 
#' @param vcf.file Path to the vcf file containing SNVs
#' @param df A bed file-like dataframe with the columns: chrom, pos, ref, alt
#' @param vcf.filter Which variants to filter for? Default is 'PASS'. Only applies if vcf.file is
#' specified
#' @param keep.chroms Which chromosomes to keep? Default is c(1:22,'X')
#' @param merge.consecutive Some vcfs report MNVs as consecutive variants. For these vcfs, such rows
#' need to be merged into one row for proper function of downstream mutSigExtractor functions.
#' @param min.muts Min number of SNVs in a PCF segment for the segment to be considered a kataegis 
#' focus
#' @param imd.thres Only PCF segments with mean intermutation distance less that this value will be
#' considered kataegis foci
#' @param cp rpart.control complexity parameter. Higher values increase overfitting of PCF
#' @param max.imd.quantile The intermutation distance at this quantile is calculated. If 
#' max.imd.quantile < imd.thres, imd.thres will be set to max.imd.quantile
#' @param output.type 'count': number of kataegis foci. 'kat_segments': simplified data of kataegis
#' segments (chrom and pos refer to the first variant in the kataegis focus). 'imd_kat': data from 
#' all kataegis associated mutations, 'imd_all': raw data. Note that the last variant from each 
#' chromosome is not reported (when output.type!='count')
#' @param verbose 
#' @param ... 
#'
#' @return An integer or dataframe
#' @export
#'
detectKataegis <- function(
   ## mutSigExtractor::variantsFromVcf) args
   vcf.file=NULL, df=NULL, vcf.filter='PASS', keep.chroms=c(1:22,'X'), 
   merge.consecutive=F,
   
   ## PCF/kataegis detection parameters
   min.muts=6, imd.thres=1000, cp=1E-6, 
   max.imd.quantile=1,
   
   output.type='count', verbose=T, ...
){
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/181206_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/171002_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/180808_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz' ## Hypermutator n_kat=1600
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/170311_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz' ## Medium hypermutator k_kat=200
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104-update3//somatics/160604_HMFreg0056_FR10302045_FR10302046_XXXXXXXX/purple/XXXXXXXX.purple.somatic.vcf.gz' ## Job fail
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor//DO217817-from-jar//purple25/DO217817T.purple.somatic.vcf.gz' ## Job fail
   # vcf.filter='PASS'
   # keep.chroms=c(1:22,'X')
   # merge.consecutive=F
   # min.muts=6
   # imd.thres=1000
   # cp=1E-5
   # 
   # verbose=T
   
   if(!is.null(vcf.file)){
      df <- vcfAsDataframe(
         vcf.file, vcf.filter='PASS', keep.chroms=keep.chroms, 
         merge.consecutive=merge.consecutive, verbose=verbose
      )
      
   } else if(!is.null(df)){
      if(c('chrom','pos','ref','alt') %in% colnames(df)){
         stop("`df` must have the colnames: chrom, pos, ref, alt")
      }
   } else {
      stop("Input must be specified to either `vcf.file` or `df`")
   }
   
   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]
   
   if(verbose){ message('Subsetting for SNVs...') }
   df <- df[nchar(df$ref)==1 & nchar(df$alt)==1,]

   if(nrow(df)==0){
      if(output.type=='count'){
         return(0)
      } else {
         return(data.frame())
      }
   }
   
   if(verbose){ message('Performing PCF on inter-mutation distances...') }
   df$chrom <- factor(df$chrom, unique(df$chrom))
   df <- df[order(df$chrom, df$pos),]
   df_split <- split(df[,c('chrom','pos')], df$chrom)
   
   imd_pcf <- do.call(rbind, lapply(df_split, function(i){
      #i=df_split[[3]]
      if(nrow(i)>=2){
         out <- data.frame(
            pos=i$pos[-nrow(i)],
            index=1:(nrow(i)-1),
            imd=i$pos[-1] - i$pos[-nrow(i)]
         )
         
         tree <- rpart::rpart(
            imd~index, data=out, 
            method='anova',
            control=rpart::rpart.control(minsplit=min.muts, cp=cp)
         )
         out$segment_id <- rle2Clusters(tree$where)
         
         ## Plotting
         # out$segment_y <- predict(tree, data.frame(index=out$index))
         # ggplot(out) +
         #    geom_point(aes(x=index, y=log10(imd))) +
         #    geom_line(aes(x=index, y=log10(segment_y)), color='red')
            
         
      } else {
         out <- data.frame(
            pos=0,
            index=1,
            imd=NA,
            segment_id=0 #,
            #segment_y=0
         )
      }
      
      out <- cbind(chrom=i[1,'chrom'], out)
      
      return(out)
   }))
   rownames(imd_pcf) <- NULL
   
   if(verbose){ message('Assigning segment numbers...') }
   ## Ensure that segments on different chromosomes have a different segment number
   imd_pcf$segment_id <- rle2Clusters(
      as.integer(paste0(as.integer(imd_pcf$chrom), imd_pcf$segment_id)) 
   )
   
   if(verbose){ message('Identifying kataegis segments...') }
   segments <- (function(){
      l <- split(
         imd_pcf[,c('segment_id','imd')], 
         #imd_pcf,
         factor(imd_pcf$segment_id, unique(imd_pcf$segment_id))
      )
      
      do.call(rbind, lapply(l, function(i){
         #i=l[[1]]
         data.frame(
            segment_id=i[1,'segment_id'],
            n_muts=nrow(i),
            mean_imd=mean(i$imd)
         )
      }))
   })()
   
   segments <- cbind(
      imd_pcf[
         match(segments$segment_id, imd_pcf$segment_id),
         c('chrom','pos')
      ],
      segments
   )
   
   # ##
   # gamma <- log(2)/median(imd_pcf$imd)
   # imd_pcf$prob <- 1 - exp(1)^(-gamma*imd_pcf$imd)
   # segments$prob_thres <- -log( 1 - (0.01/nrow(df)) ^ (1/(segments$n_muts-1)) ) / gamma
   # segments[order(segments$prob, decreasing=T),]
   # ##
   # 
   # ##
   # gamma <- log(2)/median(imd_pcf$imd)
   # segments$prob_imd <- 1 - exp(1)^(-gamma*segments$mean_imd)
   # 
   # #segments$prob_imd_consec <- nrow(df) * (1-segments$prob_imd) * segments$prob_imd^(segments$n_muts-1)
   # #segments$prob_imd_consec <- (1-segments$prob_imd) * segments$prob_imd^(segments$n_muts-1)
   # 
   # ##
   
   imd_cutoff <- min(
      imd.thres,
      quantile(imd_pcf$imd, max.imd.quantile, na.rm=T) ## To deal with hypermutators
   )
   
   kat_segments <- segments[
      with(segments,{ n_muts >= min.muts & mean_imd <= imd_cutoff })
   ,]
   imd_pcf$is_kat_locus <- imd_pcf$segment_id %in% kat_segments$segment_id
   
   if(verbose & nrow(kat_segments)==0){ message('No kataegis events found') }
   
   if(output.type=='count'){
      return(c(foci=nrow(kat_segments)))
   } else if(output.type=='kat_segments'){
      return(kat_segments)
   } else if(output.type=='imd_kat'){
      return(subset(imd_pcf, is_kat_locus))
   } else if(output.type=='imd_all'){
      return(imd_pcf)
   } else {
      stop("`output.type` must be one of the following: 'count', 'kat_segments', 'imd_kat', 'imd_all'")
   }
}

####################################################################################################
#' Make kataegis contexts
#'
#' @param vcf.snv Path to the vcf file containing SNVs
#' @param vcf.sv Path to the GRIDSS vcf file containing
#' @param vcf.filter Which variants to filter for? Default is 'PASS'. Only applies if vcf.file is
#' specified
#' @param keep.chroms Which chromosomes to keep? Default is c(1:22,'X')
#' @param merge.consecutive Some vcfs report MNVs as consecutive variants. For these vcfs, such rows
#' @param kat.foci.flank.length Number of bp on either side of kataegis to scan for SVs
#' @param verbose Show progress messages?
#'
#' @return An integer vector
#' @export
#'
extractContextsKataegis <- function(
   vcf.snv, vcf.sv, vcf.filter='PASS', keep.chroms=c(1:22,'X'), merge.consecutive=F, 
   kat.foci.flank.length=1000, verbose=T
){
   # vcf.snv='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/181206_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # vcf.sv='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/181206_HMFregXXXXXXXX/XXXXXXXX.purple.sv.ann.vcf.gz'
   
   # vcf.snv='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/161201_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # vcf.sv='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/161201_HMFregXXXXXXXX/XXXXXXXX.purple.sv.ann.vcf.gz'
   
   counts <- structure(
      rep(0,2),
      names=c('kataegis.sv_pos','kataegis.sv_neg')
   )
   
   if(verbose){ message('## Detecting kataegis') }
   kat_snv <- detectKataegis(
      vcf.snv, output.type='kat_segments', 
      vcf.filter=vcf.filter, keep.chroms=keep.chroms,
      merge.consecutive=merge.consecutive,
      verbose=verbose
   )
   
   if(nrow(kat_snv)==0){ 
      if(verbose){ message('No kataegis events found. Returning 0s') }
      return(counts) 
   }
   
   kat_snv$pos_linear <- linearizeChromPos(kat_snv$chrom, kat_snv$pos)
   
   if(verbose){ message('\n## Loading SVs') }
   variants_sv <- vcfAsDataframe(
      vcf.sv, 
      #vcf.fields=c('CHROM','POS','REF','ALT','FILTER','ID','INFO'),
      vcf.filter=vcf.filter, keep.chroms=keep.chroms, 
      verbose=verbose
   )
   
   if(verbose){ message('\n## Making kataegis contexts') }
   if(verbose){ message('Identifying unique SV breakends...') }
   if(nrow(variants_sv)!=0){
      uniq_breakends <- (function(){
         l <- regmatches(variants_sv$alt, gregexpr('\\d|\\w+:\\d+', variants_sv$alt))
         df <- as.data.frame(do.call(rbind, lapply(l, function(i){
            if(length(i) == 0){ c(NA,NA) }
            else { unlist(strsplit(i, ':')) }
         })))
         colnames(df) <- c('chrom','pos')
         df$chrom <- paste0('chr',df$chrom)
         
         df <- rbind(
            variants_sv[,c('chrom','pos')],
            df
         )
         
         df <- na.exclude(df)
         df <- unique(df)
         df$pos <- as.integer(df$pos)
         #df$chrom <- factor(df$chrom, levels=unique(variants_sv$chrom))
         #df <- df[order(df$chrom, df$pos),]
         
         return(df)
      })()
      
      uniq_breakends$pos_linear <- linearizeChromPos(uniq_breakends$chrom, uniq_breakends$pos)
      uniq_breakends <- uniq_breakends[order(uniq_breakends$pos_linear),]
      uniq_breakends <- uniq_breakends[uniq_breakends$chrom %in% kat_snv$chrom,]
   } else {
      uniq_breakends <- data.frame()
   }
   
   if(verbose){ message('Detecting nearby breakends to kataegis foci...') }
   has_nearby_breakend <- (function(){
      l <- split(kat_snv, kat_snv$segment_id)
      sapply(l, function(i){
         #i=l[[1]]
         if(nrow(uniq_breakends)!=0){
            boundary_l <- min(i$pos_linear) - kat.foci.flank.length
            boundary_r <- max(i$pos_linear) + kat.foci.flank.length
            any(boundary_l <= uniq_breakends$pos_linear & uniq_breakends$pos_linear <= boundary_r)
         } else {
            FALSE
         }
      })
   })()
   
   counts['kataegis.sv_pos'] <- sum(has_nearby_breakend)
   counts['kataegis.sv_neg'] <- sum(!has_nearby_breakend)
   
   return(counts)
}

####################################################################################################
plotRainfall <- function(imd){
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/181206_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/180808_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz' ## Hyper mutator
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/180224_HMFregXXXXXXXX/XXXXXXXX.purple.somatic.vcf.gz'
   # imd <- detectKataegis(
   #    vcf.file,
   #    output.type=2
   # )
   
   imd$chrom <- gsub('chr','',imd$chrom)
   imd$chrom <- factor(imd$chrom, unique(imd$chrom))
   #imd$is_kat_locus <- factor(imd$is_kat_locus, levels=c('FALSE','TRUE'))
   
   require(ggplot2)
   ggplot(imd, aes(x=pos, y=log10(imd))) +
      facet_grid(~chrom, space='free_x', scales='free_x') +
      geom_point(aes(color=is_kat_locus), size=0.5) +
      scale_color_manual(values=c('TRUE'='black','FALSE'='grey'), guide=F) +
      theme_bw() +
      theme(
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.spacing=unit(0, "lines"),
         panel.grid.major.x=element_blank(),
         panel.grid.minor.x=element_blank()
      )
}













