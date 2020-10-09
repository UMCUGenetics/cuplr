#' Make gene driver features from gene diplotypes table
#'
#' @param gene.diplotypes.path Path to gene_diplotypes.txt.gz (output from geneDriverAnnotator)
#' @param gene.diplotypes Alternatively the gene diplotypes dataframe can be provided as input
#' @param gene.amp.whitelist.path Path to txt file with list of selected amplified genes 
#' @param gene.def.whitelist.path Path to txt file with list of selected Deficient genes 
#'
#' @return A list: list(amp=character(), def=character())
#' @export
#'
getDrivers <- function(
   gene.diplotypes.path=NULL, gene.diplotypes=NULL,
   gene.amp.whitelist.path=GENE_AMP_WHITELIST,
   gene.def.whitelist.path=GENE_DEF_WHITELIST
){
   #dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_DR104/gene_ann_2/CPCT02300021T'
   #gene.diplotypes <- read.delim(paste0(dir,'/gene_diplotypes.txt.gz'))

   if(!is.null(gene.diplotypes.path)){
      gene.diplotypes <- read.delim(gene.diplotypes.path)
   }
   
   driver_summary <- gene.diplotypes[,c('hgnc_symbol','biall_type','hit_score','gain_type','gain_ratio_focal','gain_ratio_genome')]
   
   gene_amp_whitelist <- read.table(gene.amp.whitelist.path, header=F, stringsAsFactors=F)[,1]
   gene_def_whitelist <- read.table(gene.def.whitelist.path, header=F, stringsAsFactors=F)[,1]
   
   ## Amp
   driver_summary$gain_type[driver_summary$gain_ratio_genome < 3] <- 'none' ## select confident amp
   driver_summary$gain_type[driver_summary$gain_type!='focal'] <- 'none'
   driver_summary$gain_ratio_focal[gene.diplotypes$gain_type!='focal'] <- 0
   
   ## Def
   driver_summary$biall_type[driver_summary$hit_score < 9] <- 'none' ## select confident biall loss
   
   out <- list()
   
   out$amp <- (function(){
      df <- driver_summary[driver_summary$hgnc_symbol %in% gene_amp_whitelist,]
      #structure(df$gain_type=='focal', names=df$hgnc_symbol)
      structure(df$gain_ratio_focal, names=df$hgnc_symbol)
   })()
   
   out$def <- (function(){
      df <- driver_summary[driver_summary$hgnc_symbol %in% gene_def_whitelist,]
      structure(df$biall_type, names=df$hgnc_symbol)
   })()
   
   return(out)
}
