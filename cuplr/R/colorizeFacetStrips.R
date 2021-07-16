#' Add colors to facet strips of a ggplot object
#'
#' @param p A ggplot
#' @param colors A character vector of colors
#' @param return.gtable If TRUE will return a gtable (use this when passing the plot to cowplot to
#' preserve plot alignments). If FALSE, the gtable will be converted to a ggplot using
#' `ggplotify::as.ggplot()`
#'
#' @return A gtable or ggplot object
#' @export
#'
colorizeFacetStrips <- function(p, colors=NULL, return.gtable=T){

   if(is.null(colors)){ stop('A character vector must be provided to `colors`') } ## colors is also a base R function

   ## Convert ggplot object into raw data
   g <- ggplot_gtable(ggplot_build(p))
   names(g$grobs) <- g$layout$name

   ## Get strip ids
   strip_ids <- grep('strip-', g$layout$name, value=T)

   ## Strip ids are in the form "strip-t-1-8"
   ## Order strip ids by the last number (by default, the order is reversed)
   strip_ids_last_num <- sapply(
      strsplit(strip_ids,'-'),
      function(i){ i[length(i)] }
   )
   strip_ids <- strip_ids[order(strip_ids_last_num)]

   ## The number of strip ids is always a rectangular number (e.g. 8x7=56)
   ## If there are e.g. 54 facets, the last 2 will be empty.
   ## Check which facets are empty and exclude them
   is_empty_facet <- sapply(strip_ids, function(i){
      inherits(g$grobs[[i]],'zeroGrob')
   })
   strip_ids <- strip_ids[ !is_empty_facet ]

   if(length(colors) > length(strip_ids)){
      warning('No. of colors (',length(colors),') does not match the no. of facets (',length(strip_ids),')')
   } else if(length(colors)<length(strip_ids)){
      warning('No. of colors (',length(colors),') is less than the no. of facets (',length(strip_ids),'). Colors will be recycled')
      #colors <- c('grey','white')
      color_indexes <- rep(
         1:length(colors),
         ceiling(length(strip_ids)/length(colors))
      )
      colors <- colors[color_indexes]
   }

   ## Apply color to each strip id
   counter <- 0
   for(i in strip_ids){
      #i=strip_ids[[1]]
      counter <- counter + 1
      #message(counter,' ',i)
      object <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[object]]$gp$fill <- colors[counter]
   }

   if(return.gtable){ return(g) }

   ## Return a ggplot object
   ## However, using ggplotify messes up alignment with cowplot::plot_grid()
   ggplotify::as.ggplot(g)
}


