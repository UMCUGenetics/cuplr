#' Determine the most pathogenic diplotype effect per gene
#'
#' @param gene.diplotypes Table containing the mutation info of alleles 1 and 2
#' @param colname Which column name to take the max of?
#'
#' @return The subsetted gene.diplotypes table with a column indicating how many rows there were
#' originally per gene
#' @export
getGeneDiplotypeMaxEff <- function(gene.diplotypes){
   #gene.diplotypes=gene_diplotypes
   
   df <- gene.diplotypes
   
   ## Force: cnv_cnv as highest priority, germ+som as lowest priority
   df$diplotype_origin <- factor(
      df$diplotype_origin, 
      #c('germ_som','som_som','cnv_germ','cnv_som','cnv_cnv') 
      c('cnv_cnv','cnv_som','cnv_germ','som_som','germ_som') 
   )
   
   df <- df[
      order(
         df$ensembl_gene_id,
         -df$hit_score,
         -as.integer(df$diplotype_origin)
      )
   ,]
   df <- df[!duplicated(df$ensembl_gene_id),]
   #subset(df,hgnc_symbol=='BRCA2')
   
   ## Per gene, count number of rows with same hit_score
   
   tab <- (function(){
      l <- split(gene.diplotypes$hit_score, gene.diplotypes$ensembl_gene_id)
      unlist(lapply(l, function(i){
         if(length(i)==0){ return(0) }
         sum(i==max(i)) 
      }))
   })()
   
   df$n_max_score_biall_states <- unname(tab)[ match(df$ensembl_gene_id, names(tab)) ]
   
   df[order(df$hgnc_symbol),]
}

####################################################################################################
#' Determine biallelic hit type v2
#'
#' @param diplotypes Diplotypes dataframe
#'
#' @return A character vector of hit types
#' @export
#'
detBiallStatus <- function(gene.diplotypes){
   with(gene.diplotypes,{
      unlist(Map(function(a1.eff, a2.eff, a1.max_score, a2.max_score){
         if(a1.eff=='deep_deletion'){ 
            return('deep_deletion') 
         }
         
         if(grepl('^loh',a1.eff)){
            if(a2.max_score==5){ return('loh,mut_pathogenic') }
            if(a2.max_score==4){ return('loh,mut_likely_pathogenic') }
            if(a2.max_score==3){ return('loh,mut_vus') }
            return('loh_only')
         }
         
         if(a1.eff!='loh' & a1.max_score==5 & a2.max_score==5){ 
            return('2x_mut_pathogenic') 
         } ## germ+som and som+som cases
         
         if(a1.max_score==5 | a2.max_score==5){ return('mut_pathogenic') }
         if(a1.max_score==4 | a2.max_score==4){ return('mut_likely_pathogenic') }
         
         return('none')
      }, a1.eff, a2.eff, a1.max_score, a2.max_score, USE.NAMES=F))
   })
}

detBiallType <- function(gene.diplotypes){

   with(gene.diplotypes,{
      unlist(Map(function(a1.eff, a2.eff, a1.max_score, a2.max_score, diplotype_origin){
         
         if(a1.eff=='deep_deletion'){ return('deep_deletion') }
         
         ## LOH
         if(a1.eff=='loh' & a2.max_score>=4){ return('loh,mut') }
         if(a1.eff=='loh_arm' & a2.max_score>=4){ return('loh_arm,mut') }
         if(a1.eff=='loh_chrom' & a2.max_score>=4){ return('loh_chrom,mut') }
         
         ## germ+som and som+som cases
         if(diplotype_origin=='som_som'){
            if(a1.max_score==5 & a2.max_score==5){ return('mut,mut') }
            if(a1.max_score>=5 | a2.max_score>=5){ return('mut') }
         }
         
         if(diplotype_origin=='germ_som'){
            if(a1.max_score==5 & a2.max_score==5){ return('mut,mut') }
            if(a2.max_score>=5){ return('mut') }
         }
         
         return('none')
      }, a1.eff, a2.eff, a1.max_score, a2.max_score, diplotype_origin, USE.NAMES=F))
   })
}
