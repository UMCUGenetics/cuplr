## Load default paths on package load ================================
.onLoad <- function(libname, pkgname){

   ## Paths within the package
   pkg_constant_paths <- list(
      c('GENE_FUSION_WHITELIST','/gene_fusion_whitelist.txt'),
      
      c('GENE_DRIVER_WHITELIST','/gene_driver_whitelist.txt')
   )

   for(i in pkg_constant_paths){
      assign(
         i[1], system.file(i[2], package='featureExtractor'),
         envir=parent.env(environment())
      )
   }
}
