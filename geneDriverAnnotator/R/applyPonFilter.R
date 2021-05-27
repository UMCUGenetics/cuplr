#' Set variants to zero impact
#' 
#' @param diplotypes Diplotypes dataframe
#' @param variant.blacklist A character vector with each value in the form of {hgnc_symbol}_{hgvs_c}.
#' For example: NEIL3_c.1461-659G>A
#' @param mode Can be 'germ' or 'som'. Apply filtering to germline or somatic variants?
#' 
#' @description This function will apply the follow edits to the corresponding columns that match
#' the variants in `variant.blacklist`:
#' a*.eff -> 'none'
#' a*.max_score -> 0
#' a*.max_score_origin -> 0
#' ...where a* refers to a1 or a2, corresponding to the allele 1 or 2 columns
#'
#' @return The diplotypes dataframe with edited values
#' @export
#'
applyPonFilter <- function(diplotypes, variant.blacklist, mode='germ'){
   
   if(F){
      diplotypes=diplotypes_raw
      variant.blacklist=pon_germ
      mode='germ'
   }
   
   if(!(mode %in% c('germ','som'))){
      stop("`mode` must be 'germ' or 'som'")
   }
   
   ## CNV + mut --------------------------------
   origin.cnv_mut <- if(mode=='germ'){ 'cnv_germ' } else { 'cnv_som' }
   rows.cnv_mut <- with(diplotypes,{ 
      diplotype_origin==origin.cnv_mut & paste0(hgnc_symbol,'_',a2.hgvs_c) %in% variant.blacklist
   })
   
   diplotypes <- within(diplotypes,{
      a2.eff[rows.cnv_mut] <- 'none'
      a2.max_score[rows.cnv_mut] <- 0
      a2.max_score_origin[rows.cnv_mut] <- 'pon'
   })
   
   ## Mut + mut, allele 1 --------------------------------
   origin.mut_mut <- if(mode=='germ'){ 'germ_som' } else { 'som_som' }
   
   rows.mut_mut <- with(diplotypes,{ 
      diplotype_origin==origin.mut_mut & paste0(hgnc_symbol,'_',a1.hgvs_c) %in% variant.blacklist
   })
   
   diplotypes <- within(diplotypes,{
      a1.eff[rows.mut_mut] <- 'none'
      a1.max_score[rows.mut_mut] <- 0
      a1.max_score_origin[rows.mut_mut] <- 'pon'
   })
   
   ## Mut + mut, allele 2 --------------------------------
   if(mode=='som'){
      diplotypes <- within(diplotypes,{
         a2.eff[rows.mut_mut] <- 'none'
         a2.max_score[rows.mut_mut] <- 0
         a2.max_score_origin[rows.mut_mut] <- 'pon'
      })
   }
   
   ## Reassign biallelic status --------------------------------
   diplotypes$hit_score <- diplotypes$a1.max_score + diplotypes$a2.max_score
   diplotypes$biall_status <- detBiallStatus(diplotypes)
   diplotypes$biall_type <- detBiallType(diplotypes)
   
   return(diplotypes)
}