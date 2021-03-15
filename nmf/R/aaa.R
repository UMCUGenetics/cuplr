#========= Load default paths on package load =========#
.onLoad <- function(libname, pkgname){

   assign( 'MAIN_NMF_FUNCTIONS_PATH', system.file('../R/nmf.R', package='nmf'), envir=parent.env(environment()) )

}
