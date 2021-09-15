#' Determine the presence of gene fusions
#'
#' @param linx.fusions Path to the LINX fusion txt file, or a dataframe of the file
#' @param whitelist.path Path to gene fusion whitelist txt file
#' 
#' @return A named logical vector indicating which gene fusions are present in a sample
#' @export
#' 
getGeneFusions <- function(linx.fusions, whitelist.path=GENE_FUSION_WHITELIST){
   # if(F){
   #    linx.fusions='/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/DO51102-from-jar/linx14/DO51102T.linx.fusions.tsv'
   #    linx.fusions='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/featureExtractor/test/DO51126/DO51126T.linx.fusions.tsv'
   #    linx.fusions=cacheAndReadData(
   #       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/fusions/04_HMF_PCAWG_full/linx.fusions.merged.txt.gz'
   #    )
   # }
   
   ## Init --------------------------------
   if(is.character(linx.fusions)){
      fusions <- read.delim(linx.fusions, check.names=F)
   } else if(is.data.frame(linx.fusions)){
      fusions <- linx.fusions
      #rm(linx.fusions)
   } else {
      stop('`linx.fusions` must be a path or a dataframe')
   }
   
   ## Old version of linx use CamelCase for colnames
   ## New version uses camelCase <-- convert colnames to this format
   firstLower <- function(x) {
      substr(x, 1, 1) <- tolower(substr(x, 1, 1))
      return(x)
   }
   colnames(fusions) <- firstLower(colnames(fusions))
   
   ##
   gene_fusion_whitelist <- read.table(whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   ## Initialize output
   fusion_exists <- structure(
      rep(FALSE, length(gene_fusion_whitelist)),
      names=gene_fusion_whitelist
   )
   
   ## Select confident fusions
   fusions <- subset(fusions, likelihood %in% c('LOW','HIGH'))
   
   if(nrow(fusions)==0){ return(fusion_exists) }
   
   ## Rename fusions --------------------------------
   ## Rename one sided promiscuous fusions
   fusions$fusion_name <- (function(){
      v <- fusions$name
      
      is_promiscuous_5 <- fusions$reportedType=='PROMISCUOUS_5'
      v[is_promiscuous_5] <- paste0(fusions$geneStart[is_promiscuous_5],'_*')
      
      is_promiscuous_3 <- fusions$reportedType=='PROMISCUOUS_3'
      v[is_promiscuous_3] <- paste0('*_',fusions$geneEnd[is_promiscuous_3])
      
      return(v)
   })()
   
   ## Split fusions with 2 promiscuous genes into 2 entries
   if(sum(fusions$reportedType=='PROMISCUOUS_BOTH')>0){
      fusions <- (function(){
         df_promiscuous_both <- subset(fusions, reportedType=='PROMISCUOUS_BOTH')
         
         rbind(
            subset(fusions, reportedType!='PROMISCUOUS_BOTH'),
            within(df_promiscuous_both, fusion_name <- paste0(geneStart,'_*')),
            within(df_promiscuous_both, fusion_name <- paste0('*_',geneEnd))
         )
      })()
   }
   
   ## Group similar fusions under a common name --------------------------------
   grouping_regex <- list(
      ## HeadAndNeck_SC
      'MYB/MYBL1_NFIB'='^(MYB_NFIB|MYBL1_NFIB)$', 
      
      ## Lymphoid
      '@IGH_*'='^@IGH',
      
      ## Prostate
      '*_ETV1/4/5'='_ETV[145]$',
      # '*_ERG'='^(?!TMPRSS2_ERG).*_ERG',
      
      ## Thyroid
      '*_RET'='_RET$'
   )
   
   for(i in names(grouping_regex)){
      fusions$fusion_name[
         grepl(pattern=grouping_regex[[i]], x=fusions$fusion_name, perl=T)
      ] <- i
   }
   
   ## Output --------------------------------
   fusions <- subset(fusions, fusion_name %in% gene_fusion_whitelist)
   fusion_exists[ fusions$fusion_name ] <- TRUE
   
   return(fusion_exists)
}
