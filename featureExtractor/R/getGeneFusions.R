#' Determine the presence of gene fusions
#'
#' @param linx.fusion.path Path to LINX txt file containing fusion data
#' @param gene.fusion.whitelist.path Path to txt file with list of gene fusions
#'
#' @return A (0,1) vector
#' @export
#' 
getGeneFusions <- function(
   linx.fusion.path, gene.fusion.whitelist.path=GENE_FUSION_WHITELIST
){
   #linx.fusion.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/180502_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.fusion.tsv'
   #linx.fusion.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/161118_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.fusion.tsv'
   #linx.fusion.path='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/analysis/Arne/HMF_linx/sv-linx2/160725_HMFregXXXXXXXX/sv-linx/XXXXXXXX.linx.fusion.tsv'
   
   gene_fusion_whitelist <- read.table(gene.fusion.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   fusions <- read.delim(linx.fusion.path, check.names=F, stringsAsFactors=F)
   fusions$ReportedType[is.na(fusions$ReportedType)] <- ''
   
   if(nrow(fusions)==0){
      out <- structure(
         rep(FALSE, length(gene_fusion_whitelist)),
         names=gene_fusion_whitelist
      )
      
      return(out)
   }
      
   ## Identify intragenic fusions
   fusions <- (function(){
      l <- lapply(strsplit(fusions$Name,'_'),`[`,1:2)
      m <- do.call(rbind,l)
      colnames(m) <- c('partner_5','partner_3')
      cbind(fusions, m)
   })()
   fusions$partner_5 <- as.character(fusions$partner_5)
   fusions$partner_3 <- as.character(fusions$partner_3)
   fusions$in_same_gene <- fusions$partner_5==fusions$partner_3
   
   ## Rename promiscuous fusions
   fusions$new_name <- fusions$Name
   
   if(any(c('5P-Prom','3P-Prom') %in% fusions$ReportedType)){
      fusions <- (function(){
         df <- fusions[!fusions$in_same_gene,]
         
         df_split <- split(df, df$ReportedType)
         
         if(any('5P-Prom' %in% df$ReportedType)){
            df_split[['5P-Prom']]$new_name <- paste0(df_split[['5P-Prom']]$partner_5,'_*')
         }
         
         if(any('3P-Prom' %in% df$ReportedType)){
            df_split[['3P-Prom']]$new_name <- paste0(df_split[['3P-Prom']]$partner_3,'_*')
         }
         
         rbind(
            do.call(rbind, df_split),
            fusions[fusions$in_same_gene,]
         )
      })()
      rownames(fusions) <- NULL
   }
   
   out <- structure(
      gene_fusion_whitelist %in% fusions$new_name,
      names=gene_fusion_whitelist
   )
   
   return(out)
}
