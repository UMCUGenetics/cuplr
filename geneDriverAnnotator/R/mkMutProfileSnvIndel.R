####################################################################################################
#' Create annotations to the germline/somatic tables originating from HMF germline/somatic vcfs
#'
#' @description Primarily used to determine
#'
#' @param df.snv.indel A germline/somatic table originating from HMF germline/somatic vcfs
#' @param scoring A list containing the scoring for snpeff, clinvar, and enigma annotation values
#' @param keep.only.first.eff Only keep first snpeff_eff (items are separated by '&')?
#' @param include.hotspots If TRUE, identified hotspot mutations will be considered in the final
#' scoring
#' @param filter.no.impact.variants Variants with a score of 0 will be assigned 'no_impact_variant'
#' to snpeff_eff
#' @param verbose Show messages?
#'
#' @return The original input dataframe with annotation columns
#' @export
#'
mkMutProfileSnvIndel <- function(
   df.snv.indel, 
   scoring=SCORING_MUT,
   keep.only.first.eff=T,
   include.hotspots=F,
   filter.no.impact.variants=T,
   verbose=T
){
   
   #--------- Sanity checks ---------#
   if(nrow(df.snv.indel)==0){ 
      return(data.frame()) 
   }
   #df.snv.indel=input_tables$som_txt
   #df.snv.indel=input_tables$germ_txt
   
   #--------- Score annotations ---------#
   if(verbose){ message('Getting snpEff scores...') }
   snpeff_eff <- sapply(strsplit(df.snv.indel$snpeff_eff,'&'),`[[`,1)
   
   snpeff_score <- unname(scoring$snpeff)[ match(snpeff_eff,names(scoring$snpeff)) ]
   snpeff_score[is.na(snpeff_score)] <- 0 ## If annotation not found in snpeff table, return 0
   
   if(verbose){ message('Getting ClinVar scores...') }
   clinvar_score <- unname(scoring$clinvar)[ match(df.snv.indel$clinvar_sig,names(scoring$clinvar)) ]
   clinvar_score[is.na(clinvar_score)] <- 0
   
   if(verbose){ message('Assigning hotspot score...') }
   if(include.hotspots){
      hotspot_score <- ifelse(df.snv.indel$is_hotspot_mut,5,0)
   } else {
      hotspot_score <- rep(0, nrow(df.snv.indel))
   }
   
   #--------- Calculate max score and which database it came from ---------#
   if(verbose){ message('Calculating max scores...') }
   sig_scores <- data.frame(clinvar_score, hotspot_score, snpeff_score, stringsAsFactors=F)
   
   sig_scores$max_score <- unlist(Map(function(clinvar_score, hotspot_score, snpeff_score){
      if(clinvar_score != 0){ return(clinvar_score) }
      if(hotspot_score != 0){ return(hotspot_score) }
      if(snpeff_score != 0){ return(snpeff_score) }
      return(0)
   }, clinvar_score, hotspot_score, snpeff_score, USE.NAMES=F))
   
   sig_scores$max_score_origin <- unlist(Map(function(clinvar_score, hotspot_score, snpeff_score){
      if(clinvar_score != 0){ return('clinvar') }
      if(hotspot_score != 0){ return('hotspot') }
      if(snpeff_score != 0){ return('snpeff') }
      return('none')
   }, clinvar_score, hotspot_score, snpeff_score, USE.NAMES=F))
   
   #--------- Output ---------#
   out <- cbind(df.snv.indel, sig_scores)
   
   ## Deal with multiple sneff effs
   if(keep.only.first.eff){
      if(verbose){ message('Keeping only first snpeff_eff...') }
      out$snpeff_eff <- gsub('&.+$','',out$snpeff_eff)
   }
   
   if(filter.no.impact.variants){
      out$snpeff_eff[ 
         !(out$snpeff_eff %in% names(scoring$snpeff)) & out$max_score==0
      ] <- 'no_impact_variant'
   }
   
   return(out)
}
