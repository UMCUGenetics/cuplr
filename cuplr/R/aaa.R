## Load default paths on package load ================================

.onLoad <- function(libname, pkgname){

   ## Paths within the package
   pkg_constant_paths <- list(
      c('MODEL_PATH','/model.rds'),
      c('PROB_CAL_CURVES_PATH','/prob_calib_curves.txt')
   )

   for(i in pkg_constant_paths){
      assign(
         i[1], system.file(i[2], package='cuplr'),
         envir=parent.env(environment())
      )
   }
}

# if(F){
#    model_dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.17a_noClonalSplit/'
#    pkg_dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/'
#
#    file.copy(
#       paste0(model_dir,'/final/model.rds'),
#       paste0(pkg_dir,'/cuplr/inst/model.rds'),
#       overwrite=T
#    )
#
#    file.copy(
#       paste0(model_dir,'/report/prob_calib_curves.txt'),
#       paste0(pkg_dir,'/cuplr/inst/prob_calib_curves.txt'),
#       overwrite=T
#    )
# }

