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
   vcf.file=NULL, df=NULL, sample.name=NULL, vcf.filter='PASS', keep.chroms=c(1:22,'X'),
   ref.genome=mutSigExtractor::DEFAULT_GENOME, #clonal.variants.only=T,
   
   ## PCF/kataegis detection parameters
   min.muts=6, imd.thres=1000, cp=1E-6, 
   max.imd.quantile=1,
   
   output.type='count', verbose=T, ...
){

   if(!is.null(vcf.file)){
      
      df <- mutSigExtractor::variantsFromVcf(
         vcf.file, vcf.fields=c(1,2,4,5,7,8),
         vcf.filter=vcf.filter, keep.chroms=keep.chroms,
         ref.genome=ref.genome, verbose=verbose,
         ...
      )
   }
   if(is.null(df)){ stop("Input must be specified to either `vcf.file` or `df`") }
   
   if(nrow(df)!=0){
      df_colnames <- c('chrom','pos','ref','alt')
      if(!(identical(colnames(df)[1:4], df_colnames))){
         warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
         colnames(df)[1:4] <- df_colnames
      }
      
      if(verbose){ message('Removing rows with multiple ALT sequences...') }
      df <- df[!grepl(',',df$alt),]
      
      if(verbose){ message('Subsetting for SNVs...') }
      df <- df[nchar(df$ref)==1 & nchar(df$alt)==1,]
   }
   
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
   
   ## Helper function
   rle2Clusters <- function(x){
      rle_out <- rle(x)
      unlist(lapply(1:length(rle_out$lengths), function(i){
         rep(i, rle_out$lengths[i])
      }))
   }
   
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

