## Load default paths on package load ================================
.onLoad <- function(libname, pkgname){

   ## Paths within the package
   pkg_constant_paths <- list(
      c('GENE_FUSION_WHITELIST','/gene_fusion_whitelist.txt'),
      c('GENE_DRIVER_WHITELIST','/gene_driver_whitelist.txt'),
      c('RMD_BINS_HG19_PATH','/rmd_bins.hg19.txt.gz'),
      c('RMD_BINS_HG38_PATH','/rmd_bins.hg38.hg19_liftover.txt.gz')
   )

   for(i in pkg_constant_paths){
      assign(
         i[1], system.file(i[2], package='featureExtractor'),
         envir=parent.env(environment())
      )
   }
   
   ## Set reference genome
   setGenome('hg19', verbose=1)
}

# if(F){
#    pkg_dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/'
#    
#    file.copy(
#       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/fusions/05_newLabels/fusion_whitelist.txt',
#       paste0(pkg_dir,'/cuplr/inst/gene_fusion_whitelist.txt'),
#       overwrite=T
#    )
#    
#    file.copy(
#       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/fusions/06_newLabels/fusion_whitelist.txt',
#       paste0(pkg_dir,'/cuplr/inst/gene_driver_whitelist.txt'),
#       overwrite=T
#    )
# }