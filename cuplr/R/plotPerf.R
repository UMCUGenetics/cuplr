#' Plot a heatmap from a matrix
#'
#' @param m A numeric matrix
#' @param x.lab x-axis labels
#' @param y.lab y-axis labels
#' @param legend.name Legend title
#' @param show.labels Show raw values within each cell?
#' @param palette Name of the brewer color palette
#' @param palette.direction Color palette direction. Can be 1 (forward) or -1 (reverse).
#' @param invert.y Invert the y-axis?
#' @param x.title.position Position of the x axis title
#'
#' @return A ggplot object
#' @export
#'
plotHeatmapFromMatrix <- function(
   m, x.lab=NULL, y.lab=NULL, legend.name='value',
   show.labels=F, custom.labels=NULL, palette='YlGnBu', palette.direction=-1, invert.y=F,
   x.title.position='bottom'
){

   #m=t(imp)
   if(!is.matrix(m) & !is.numeric(m)){ stop('`m` must be a numeric matrix') }
   m_melt <- as.data.frame(reshape2::melt(m))
   colnames(m_melt)[1:2] <- c('y','x')
   m_melt$x <- factor(m_melt$x, levels=unique(m_melt$x))
   m_melt$y <- factor(m_melt$y, levels=unique(m_melt$y))
   m_melt$label <- m_melt$value

   if(invert.y){
      m_melt$y <- factor(m_melt$y, rev(levels(m_melt$y)))
   }

   p <- ggplot(m_melt, aes(x, y)) +
      scale_x_discrete(expand=c(0,0), position=x.title.position) +
      scale_y_discrete(expand=c(0,0)) +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x.bottom=element_text(angle=90, vjust=0.5, hjust=1),
         axis.text.x.top=element_text(angle=90, vjust=0.5, hjust=0),
         axis.ticks=element_blank()
      )

   if(palette=='none'){
      p <- p + geom_tile(fill='white', color='grey')
   } else {
      p <- p +
         geom_tile(aes(fill=value), color='grey') +
         scale_fill_distiller(
            name=legend.name,
            palette=palette, direction=palette.direction,
            guide=guide_colorbar(
               frame.colour='black', ticks.colour='black',
               direction='horizontal', title.position='top', reverse=T, barwidth=4, barheight=1,
               label.theme=element_text(angle=90, vjust=0.5, size=10)
            )
         )
   }

   if(!is.null(custom.labels)){
      if(is.matrix(custom.labels) | is.data.frame(custom.labels)){
         if(!identical(dim(custom.labels),dim(m))){
            stop('custom.labels when a matrix must have the same dims as m: ',nrow(m),'x',ncol(m))
         }
         m_melt$label <- melt(as.matrix(custom.labels))$value
      } else {
         if(length(custom.labels) != nrow(m_melt)){
            stop('custom.labels when a vector must be the same length as the number of observations: ',nrow(m_melt))
         }
         m_melt$label <- custom.labels
      }
   }

   if(show.labels){
      p <- p + geom_text(data=m_melt, aes(label=label), size=3.5)
   }

   if(is.null(x.lab)){
      p <- p + theme(axis.title.x.top=element_blank(), axis.title.x.bottom=element_blank())
   } else {
      p <- p + xlab(x.lab)
   }

   if(is.null(y.lab)){
      p <- p + theme(axis.title.y=element_blank())
   } else {
      p <- p + ylab(y.lab)
   }

   return(p)
}

####################################################################################################
#' Plot heatmap of classifier performance
#'
#' @description Creates two plots. Upper plot shows performance metrics per cancer type and lower
#' plot shows the number or % of samples correctly/incorrectly classified
#'
#' @param actual A vector of the actual classes
#' @param predicted A vector of the predicted classes
#' @param cv_out Alternative input to `actual` and `predicted`. Cross validation output in the
#' form of a list that has 'imp' object (e.g. cv_out[[1]]$test_set$actual and
#' cv_out[[1]]$test_set$predicted)
#' @param rel.heights A numeric vector of length 3. Relative heights of the plots
#' @param sort Sort cancer type by number of % sample of correctly classified
#' @param rel.values In lower plot, show absolute number or %
#' @param metrics A character vector indicating which performance metrics to show in the upper plot.
#' See documentation for mltoolkit to see which metrics are available
#' @param show.weighted.mean In upper plot, calculated the weighted mean or normal mean?
#'
#' @return A cowplot object
#' @export
#'
plotPerfHeatmap <- function(
   actual=NULL, predicted=NULL, cv_out=NULL,
   rel.heights=c(perf=0.3, counts=0.15, heatmap=1),

   ## Confusion heatmap
   sort=F, rel.values=T,

   ## Performance metrics
   metrics=c('f1','prec','tpr'),
   show.weighted.mean=T
){
   if(F){
      test_set=readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.13a_HMF_PCAWG/test_set.rds')
      actual=test_set$actual
      predicted=test_set$predicted
   }

   require(ggplot2)
   if(!is.null(cv_out)){
      actual <- unlist(lapply(cv_out,function(i){ i$test_set$actual }))
      predicted <- unlist(lapply(cv_out,function(i){ i$test_set$predicted }))
   }

   ## Make levels the same between actual and predicted
   classes <- sort(unique(c(
      as.character(predicted),
      as.character(actual)
   )))

   classes <- c(
      classes[classes %in% predicted],
      classes[!(classes %in% predicted)]
   )

   actual <- factor(actual, classes)
   predicted <- factor(predicted, classes)

   ## Large confusion matrix ----------------------------------------------------------------
   tab <- table(predicted, actual)

   if(rel.values){
      tab <- apply(tab,2,function(i){ i/sum(i) })
      tab <- round(tab,2)
   }

   class_order <- if(sort){ names(sort(diag(tab), decreasing=T)) } else { names(diag(tab)) }
   tab <- tab[class_order,class_order]

   tab <- cbind(NA,tab) ## Add dummy column for mean performance
   tab <- tab[levels(droplevels(predicted)),] ## remove non-existing prediction levels

   plots <- list()

   plots$heatmap <- plotHeatmapFromMatrix(
      tab, show.labels=T, x.lab='Actual\n(Columns: prop. misclassified as which class)', y.lab='Predicted', invert.y=T,
      #palette=if(rel.values){ 'YlGnBu' } else { 'none' },
      legend.name=if(rel.values){ 'Column fraction' } else { 'Counts' }
   ) +
      theme(
         axis.text.x.bottom=element_blank()
      )

   ## Samples per class ----------------------------------------------------------------
   class_counts <- table(actual)
   class_counts <- class_counts[class_order]
   class_counts <- c(total=length(actual), class_counts)

   correct_counts <- unlist(lapply(split(data.frame(actual, predicted), actual), function(i){
      sum(i$actual==i$predicted)
   }))
   correct_counts <- correct_counts[class_order]
   correct_counts <- c(total=sum(actual==predicted),correct_counts)

   counts <- data.frame(total=class_counts, correct=correct_counts)
   counts$incorrect <- counts$total - counts$correct
   #counts <- counts[,ncol(counts):1]
   #counts <- forceDfOrder(counts)

   plots$counts <-
      plotHeatmapFromMatrix(t(as.matrix(counts)), palette='none', show.labels=T, x.title.position='top') +
      geom_hline(yintercept=1.5, size=0.5) +
      theme(
         axis.text.x.bottom=element_blank(),
         axis.title=element_blank()
      )

   ## Performance stats ----------------------------------------------------------------
   confusion <- mltoolkit::confusionMatrix(
      predicted=mltoolkit::oneHotEncode(predicted),
      actual=actual,
      simplify=T
   )

   perf <- mltoolkit::calcPerf(confusion, metrics)
   rownames(perf) <- perf[,1]; perf[,1] <- NULL

   perf <- perf[class_order,,drop=F]

   perf_ss <- perf[rownames(perf) %in% predicted,,drop=F]

   ## Summary stats
   if(show.weighted.mean){
      perf_summary <- rbind(
         OVERALL=apply(perf_ss,2,function(i){
            weights <- table(actual)/length(actual)
            weights <- weights[names(i)]
            weighted.mean(i, weights)
         })
      )
   } else {
      perf_summary <- rbind(
         OVERALL = apply(perf_ss,2,mean)
      )
   }
   perf <- rbind(perf_summary, perf)
   perf <- round(perf, 2)

   plots$perf <- plotHeatmapFromMatrix(
         t(perf), palette='RdYlGn', palette.direction=1, show.labels=T,
         y.lab='Perf.', legend.name='Metric value', x.title.position='top'
      ) +
      guides(fill=F) +
      theme(
         axis.text.x.bottom=element_blank(),
         axis.title=element_blank()
      )

   ## Combine --------------------------------
   plots <- rev(plots)

   sel_plots <- names(rel.heights!=0)
   plots <- plots[sel_plots]
   rel_heights <- rel.heights[sel_plots]

   plots <- lapply(plots, function(i){
      i + theme(
         axis.text.x.bottom=element_blank(),
         axis.text.x.top=element_blank()
      )
   })

   plots[[1]] <- plots[[1]] +
      theme(
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5)
      )

   plots[[length(plots)]] <- plots[[length(plots)]] +
      theme(
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x.bottom=element_text()
      )

   cowplot::plot_grid(
      plotlist=plots,
      ncol=1, align='v', axis='tblr', rel_heights=rel_heights
   )

}


####################################################################################################
#' Plot false negative rate
#'
#' @description The false negative rate can also be interpreted as the cumulative fraction of
#' samples per class that have a given probability.
#'
#' @param actual A vector of the actual classes
#' @param prob A matrix of probabilities, where rows are the samples and columns are the cancer
#' types
#'
#' @return A ggplot object
#' @export
#'
plotFnr <- function(actual, prob){

   confusion <- mltoolkit::confusionMatrix(
      predicted=prob,
      actual=actual,
      cutoff.interval=0.1
   )

   perf <- do.call(rbind, lapply(names(confusion), function(i){
      df <- confusion[[i]]
      df <- cbind(
         df, mltoolkit::calcPerf(df,metrics=c('tpr','fnr','prec','f1'), add.start.end.values=F)
      )
      df <- cbind(class=i,df)
      return(df)
   }))

   m <- reshape2::dcast( perf[,c('class','cutoff','fn')], class~cutoff, value.var='fn' )
   rownames(m) <- m[,1]; m[,1] <- NULL
   m <- as.matrix(m)

   m_rel <- round(m/m[,ncol(m)],2)

   custom_labels <- paste0(
      format(round(reshape2::melt(m_rel)$value, 2), nsmall=2),
      ' (',reshape2::melt(m)$value,')'
   )

   plotHeatmapFromMatrix(
      m_rel, invert.y=T, x.lab='Probability', y.lab='Cancer type',
      show.labels=T, custom.labels=custom_labels, legend.name='Fraction'
   ) +
      #ggtitle('Per cancer type, fraction/counts of samples <=probability') +
      ggtitle(expression("Per cancer type, fraction/counts of samples with prob."<="indicated prob.")) +
      theme(
         axis.text.x.bottom=element_text(angle=0, hjust=0.5)
      )

}
