#' Plot regional mutational density profiles
#'
#' @param m A matrix where columns are RMD bins and rows are samples/signature profiles
#' @param top.n Show labels for the features with the highest peaks
#' @param plot.title Plot title
#' @param y.lab y axis title
#'
#' @return A ggplot object
#' @export
#'
plotRmdProfile <- function(m, top.n=15, plot.title='RMD profile', y.lab='Probability'){
   #m=model$rmd_sig_profiles
   #m <- t(m)

   require(ggplot2)
   require(ggrepel)

   #if(!is.matrix(m)){ stop('`m` must be a matrix') }

   if(is.vector(m)){
      m <- matrix(m, nrow=1, dimnames=list(plot.title,names(m)))
   }

   ## Init plot data
   pd <- reshape2::melt(as.matrix(m))
   colnames(pd) <- c('sample','feature','probability')

   pd$sample <- factor(pd$sample, unique(pd$sample))

   pd$chrom <- regmatches(as.character(pd$feature), regexpr('((\\d+)|[XYM])[pq]+', as.character(pd$feature)))
   pd$chrom <- gsub('p|q','',pd$chrom)
   pd$chrom <- factor(pd$chrom, naturalsort::naturalsort(unique(pd$chrom)))

   pd$start <- as.numeric(sapply(strsplit(as.character(pd$feature),'_'), `[[`, 2))

   ## Calculate linear chrom position
   uniq_bins <- unique(pd[,c('chrom','start')])
   chrom_lengths <- sapply(split(uniq_bins, uniq_bins$chrom), function(i){
      max(i$start) + 1
   })
   #chrom_lengths <- commonUtils::CHROM_LENGTHS_HG19
   names(chrom_lengths) <- gsub('chr','',names(chrom_lengths))
   #chrom_lengths <- chrom_lengths/1e6

   offsets <- cumsum(as.numeric(chrom_lengths))
   offsets <- c(0,offsets[-length(offsets)])
   names(offsets) <- names(chrom_lengths)

   pd$start_linear <- pd$start + offsets[ as.character(pd$chrom) ]

   ## Assign labels to top bins
   if(!is.null(top.n)){
      if(top.n>30){ stop('top.n too large') }

      top_n <- top.n
      if(top.n > ncol(m)){ top_n <- ncol(m) }

      pd_split <- split(pd, pd$sample)
      pd <- do.call(rbind, lapply(pd_split, function(i){
         #i=pd_split[[1]]
         i <- i[order(i$probability, decreasing=T),]
         i$is_top <- FALSE
         i$is_top[1:top.n] <- TRUE
         return(i)
      }))
      rownames(pd) <- NULL
   }
   pd$label <- as.character(pd$feature)
   pd$label[!pd$is_top] <- ''

   ## Plot params
   chrom_ends_linear <- offsets[levels(pd$chrom)]

   chrom_mids_linear <- chrom_lengths/2 + offsets
   chrom_mids_linear <- chrom_mids_linear[levels(pd$chrom)]

   color_pal <- rep(RColorBrewer::brewer.pal(8, 'Dark2'),3)
   color_pal <- structure(
      color_pal[1:length(levels(pd$chrom))],
      names=levels(pd$chrom)
   )

   ## Main
   l <- lapply(levels(pd$sample), function(i){
      #i='Lymphoid.1'
      pd_ss <- subset(pd, sample==i)
      ggplot(pd_ss, aes(x=start_linear, y=probability)) +
         #facet_grid(~chrom, scales='free_x', space='free_x') +

         geom_step(aes(color=chrom),size=0.4) +
         scale_color_manual(values=color_pal, guide=F) +

         geom_text_repel(
            aes(label=label),
            #nudge_y=max(pd$probability)*0.005
            min.segment.length=0

         ) +

         geom_vline(xintercept=chrom_ends_linear, color='grey', size=0.3) +
         scale_x_continuous(
            breaks=chrom_mids_linear, labels=names(chrom_mids_linear),
            expand=c(0,0), name='chrom:pos'
         ) +

         ggtitle(i) +
         ylab(y.lab) +

         theme_bw() +
         theme(
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.spacing = unit(0, "lines")
            #axis.text.x=element_blank()
         )
   })

   names(l) <- levels(pd$sample)

   if(length(l)==1){ l[[1]] } else { l }
}

