#' Perform an enrichment (fisher test) analysis across multiple classes
#'
#' @param df An input dataframe
#' @param feature.colname The colname of the feature (e.g. 'gene')
#' @param class.colname The colname of the classes (default: 'cancer_type')
#' @param v.classes A vector of classes from the cohort. Used to determine the
#' class background frequencies
#' @param verbose Show progress messages?
#'
#' @return A dataframe containing the fisher test data
#' @export
#'
#' @examples
#' ## Load data
#' metadata <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/metadata/cancer_type_labels.txt')
#' metadata <- metadata[metadata$cohort=='HMF',]
#' driver_catalog <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/select_drivers//driver.catalog.merged.txt.gz')
#' driver_catalog$cancer_type <- metadata[match(driver_catalog$sample, metadata$sample),'cancer_type']
#' amp_catalog <- driver_catalog[driver_catalog$driver=='AMP',]
#' 
#' ##
#' enrichment <- calcCancerTypeEnrichment(
#'   amp_catalog, feature.colname='gene',
#'   v.classes=metadata$cancer_type
#' )
calcCancerTypeEnrichment <- function(
   df, feature.colname, class.colname='cancer_type',
   v.classes, verbose=F
){
   #df=amp_catalog
   #class.colname='cancer_type'
   #feature.colname='gene'
   #v.classes=metadata$cancer_type
   
   ##--------------------------------
   if(verbose){ message('Calculating occurrence of feature per cancer type...') }
   l <- split(
      df[,c(class.colname, feature.colname)], 
      df[,feature.colname]
   )
   
   pos_freqs <- do.call(rbind, lapply(l, function(i){
      #i=l[[1]]
      tab <- table(i[,class.colname])
      data.frame(
         feature=i[1,feature.colname],
         class=names(tab),
         pos_freq=as.integer(tab),
         pos_total=nrow(i)
      )
   }))
   
   ##--------------------------------
   if(verbose){ message('Calculating cancer type frequencies...') }
   background_freqs <- (function(){
      tab <- table(v.classes)
      freqs <- data.frame(
         class=names(tab),
         freq=as.integer(tab),
         total=sum(tab)
      )
      freqs$freq_rel <- freqs$freq / freqs$total
      return(freqs)
   })()
   
   neg_freqs <- background_freqs[
      match(pos_freqs$class, background_freqs$class),
      c('freq','total')
   ]
   colnames(neg_freqs) <- c('neg_freq','neg_total')
   
   out <- cbind(pos_freqs, neg_freqs)
   rownames(out) <- NULL
   
   out$pos_rel_req <- out$pos_freq/out$pos_total
   out$neg_rel_req <- out$neg_freq/out$neg_total
   
   ##--------------------------------
   if(verbose){ message('Performing fisher tests...') }
   
   conting <- as.matrix(out[,c('pos_freq','pos_total','neg_freq','neg_total')])
   colnames(conting) <- NULL
   
   if(verbose){
      counter <- 0
      pb <- txtProgressBar(max=nrow(conting), style=3)
   }
   
   out$p_value <- apply(conting, 1, function(i){
      #i=conting[1,]
      if(verbose){
         counter <<- counter + 1
         setTxtProgressBar(pb, counter)
      }
      
      m <- matrix(i, nrow=2)
      fisher.test(m, alternative='greater')$p.value
   })
   out$q_value <- p.adjust(out$p_value, method='bonferroni')
   
   out <- out[order(out$p_value),]
   out$minus_log10_pvalue <- -log10(out$p_value)
   out$minus_log10_qvalue <- -log10(out$q_value)
   
   return(out)
}

################################################################################
#' Plot enrichment analysis
#'
#' @param df The output from calcCancerTypeEnrichment()
#' @param keep.y.order Force the original order of the features?
#' @param invert.y Invert the y-axis?
#' @param show.top.n.features Show only n features (subset y axis)
#' @param show.qvalues Show q values instead of p values?
#'
#' @return A ggplot object
#' @export
#'
plotCancerTypeEnrichment <- function(
   df, keep.y.order=F, invert.y=T, show.top.n.features=NULL, show.qvalues=F
){
   #df=amp_enrichment_ss
   
   ## Subsetting --------------------------------
   if(!is.null(show.top.n.features)){
      n_features <- length(unique(df$feature))
      
      top_n <- 
         if(show.top.n.features > n_features){
            n_features
         } else {
            show.top.n.features
         }
      
      feature_whitelist <- df$feature[!duplicated(df$feature)][1:top_n]
      df <- df[df$feature %in% feature_whitelist,]
   }
   
   ## y-axis ordering  --------------------------------
   if(keep.y.order){
      df$feature <- factor(df$feature, levels=unique(df$feature))
   } else {
      df$feature <- as.factor(df$feature)
   }
   
   if(invert.y){
      df$feature <- factor(df$feature, levels=rev(levels(df$feature)))
   }
   
   ## Plotting --------------------------------
   require(ggplot2)
   p <- ggplot(df, aes(y=feature, x=class))
   
   if(!show.qvalues){
      p <- p +
         geom_point(aes(size=minus_log10_pvalue, fill=minus_log10_pvalue), shape=21) +
         scale_fill_distiller(palette='Spectral', name='-log10(pvalue)')
   } else {
      p <- p +
         geom_point(aes(size=minus_log10_qvalue, fill=minus_log10_qvalue), shape=21) +
         scale_fill_distiller(palette='Spectral', name='-log10(qvalue)')
   }

   p +
      guides(size=F) +
      ylab('Feature') +
      xlab('Cancer type') +
      theme_bw() +
      theme(
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
      )
}
