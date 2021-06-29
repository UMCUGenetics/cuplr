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

