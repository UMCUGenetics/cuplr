## Load default paths on package load ================================

# .onLoad <- function(libname, pkgname){
#
#    ## Paths within the package
#    pkg_constant_paths <- list(
#       c('CUPLR','/model.rds')
#    )
#
#    for(i in pkg_constant_paths){
#       assign(
#          i[1], system.file(i[2], package='cuplr'),
#          envir=parent.env(environment())
#       )
#    }
# }
#
# if(F){
#    file.copy(
#       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.09c_originalFeatures/model.rds',
#       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/inst/model.rds'
#    )
# }

