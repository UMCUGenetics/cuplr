#========= Load default paths on package load =========#
.onLoad <- function(libname, pkgname){

   ## Paths within the package
   pkg_constant_paths <- list(
      c('CYTOBANDS','/cytoBand.txt.gz'),
      
      c('VIRAL_HOST_REF','/viral_host_ref_ann.txt'),
      c('GENE_FUSION_WHITELIST','/gene_fusion_whitelist.txt'),
      
      c('GENE_AMP_WHITELIST','/gene_amp_whitelist.txt'),
      c('GENE_BIALL_WHITELIST','/gene_biallel_whitelist.txt'),
      c('GENE_MONOALL_WHITELIST','/gene_monoallel_whitelist.txt'),
      c('GENE_MUT_WHITELIST','/gene_mut_whitelist.txt'),
      
      c('REP_ELEM_WHITELIST','/rep_elem_whitelist.txt'),
      
      c('RESOLVED_TYPE_ANNOTATIONS','/resolved_type_annotations.txt')
   )

   for(i in pkg_constant_paths){
      assign(
         i[1], system.file(i[2], package='hmfSvFeatureExtractor'),
         envir=parent.env(environment())
      )
   }
}

if(F){
   devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/')
}