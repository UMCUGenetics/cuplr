#' Export a list of plots as a PDF with a PNG per page
#'
#' @param plot.list A list of plots
#' @param out.path Output PDF path
#' @param page.width Page width in inches
#' @param page.height Page height in inches
#' @param res Page resolution
#'
#' @export
#'
exportPlotsAsPngPdf <- function(
   plot.list, out.path,
   page.width=10, page.height=4, res=300
){

   if(!grepl('pdf$',out.path)){ stop('Output file extension must be a pdf') }

   message('Making tmp dir...')
   tmp_dir <- paste0(dirname(out.path),'/tmp_',basename(out.path),'/')
   if(dir.exists(tmp_dir)){
      system(sprintf('rm -r %s', tmp_dir))
   }
   system(sprintf('mkdir -p %s', tmp_dir))

   message('Creating PNGs...')
   for(i in 1:length(plot.list)){

      message('[',i,'/',length(plot.list),'] ',names(plot.list)[[i]])

      png_path <- paste0(
         tmp_dir,'/',
         formatC(i, width=4, format="d", flag="0"),
         '.png'
      )

      png(png_path, page.width, page.height, res=res, units='in')
      plot(plot.list[[i]])
      dev.off()
   }

   message('Merging PNGs into one PDF...')
   system(sprintf(
      'convert %s %s',
      paste0(tmp_dir,'/*.png'),
      out.path
   ))

   message('Removing tmp dir...')
   system(sprintf('rm -r %s', tmp_dir))
}
