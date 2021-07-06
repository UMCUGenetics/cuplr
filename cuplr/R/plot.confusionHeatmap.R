#' Plot a heatmap from a matrix
#'
#' @param m A numeric matrix
#' @param x.lab x-axis labels
#' @param y.lab y-axis labels
#' @param legend.name Legend title
#' @param show.labels Show value labels within each cell?
#' @param custom.labels A matrix of the same dimensions as `m`. When `show.labels`==TRUE, these
#' custom value labels will be shown instead
#' @param label.size Size of the value labels
#' @param palette Name of the brewer color palette
#' @param palette.direction Color palette direction. Can be 1 (forward) or -1 (reverse).
#' @param invert.y Invert the y-axis?
#' @param x.title.position Position of the x axis title
#'
#' @return A ggplot object
#' @export
#'
heatmapFromMatrix <- function(
   m, x.lab=NULL, y.lab=NULL, legend.name='value',
   show.labels=F, custom.labels=NULL, label.size=2.5,
   palette='YlGnBu', palette.direction=-1, fill.limits=NULL,
   invert.y=F, x.title.position='bottom'
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
            palette=palette, direction=palette.direction, limits=fill.limits,
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
      p <- p + geom_text(data=m_melt, aes(label=label), size=label.size)
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
#' @param actual A factor vector of the actual classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param predicted A factor vector of the predicted classes
#' @param which.plots Which plots to show? A character vector with one or more of the following
#' values: 'perf','counts','confusion'.
#' @param rel.heights A numeric vector of length 3. Relative heights of the plots
#' @param sort.classes Sort cancer type by number of % sample of correctly classified
#' @param rel.values In lower plot, show absolute number or %
#' @param predicted.classes.only Only show classes that are present in predicted
#' @param metrics A character vector indicating which performance metrics to show in the upper plot.
#' See documentation for mltoolkit to see which metrics are available
#' @param plot.title Plot title
#' @param label.size Size of the value labels
#'
#' @return A cowplot object
#' @export
#'
confusionHeatmap <- function(
   actual=NULL, probs=NULL, predicted=NULL,
   which.plots=c('counts','perf','confusion'),
   rel.heights=c(counts=0.3, perf=0.1, confusion=1),

   ## Confusion heatmap
   sort.classes=F, rel.values=T, predicted.classes.only=F,

   ## Performance metrics
   metrics=c('precision','f1'), show.weighted.mean=T,

   ## Plotting args
   plot.title=NULL, label.size=2.5
){
   if(F){
      report=pred_reports$CV
      actual=report$class_actual
      probs=report$prob_scaled

      sort.classes=F; rel.values=T; predicted.classes.only=F;
      which.plots=c('counts','perf','confusion')
      rel.heights=c(counts=0.3, perf=0.1, confusion=1)

      metrics=c('precision','f1','true positive rate')
      show.weighted.mean=T
      plot.title=NULL; label.size=2.5
   }

   ## Init ----------------------------
   require(ggplot2)

   if(!is.null(probs)){
      predicted <- factor(
         colnames(probs)[ max.col(probs) ],
         colnames(probs)
      )
   }

   if(!is.factor(actual)){ stop('`actual` must be a factor') }
   if(!is.factor(predicted)){ stop('`predicted` must be a factor') }

   ## Only show classes that are present in predicted
   if(!is.null(probs)){
      predicted_classes <- colnames(probs)
   } else {
      predicted_classes <- levels(predicted)
   }

   if(predicted.classes.only){
      is_predicted_class <- actual %in% predicted_classes
      actual <- actual[is_predicted_class]
      predicted <- predicted[is_predicted_class]

   } else {
      ## Make levels the same between actual and predicted
      classes <- sort(unique(c(
         predicted_classes,
         levels(actual)
      )))

      classes <- c(
         classes[classes %in% predicted_classes],
         classes[!(classes %in% levels(actual))]
      )
   }

   actual <- factor(actual, classes)
   predicted <- factor(predicted, classes)

   ## Initialize output
   plots <- list()

   ## Large confusion matrix ----------------------------
   tab <- table(predicted, actual)

   overall_acc <- sum( as.character(predicted)==as.character(actual) )
   if(rel.values){
      tab <- apply(tab,2,function(i){ i/sum(i) })
      tab <- round(tab,2)

      overall_acc <- round(overall_acc / length(actual), 2)
   }

   ## Sort classes
   class_order <- if(sort.classes){
      names(sort(diag(tab), decreasing=T))
   } else {
      names(diag(tab))
   }
   tab <- tab[class_order,class_order]

   ## Only keep predicted classes on prediction axis
   tab <- tab[predicted_classes,]

   ## Add dummy row and column for overall performance
   tab <- cbind(OVERALL=NA,tab)
   tab <- rbind(OVERALL=NA,tab)
   tab[1,1] <- overall_acc

   if('confusion' %in% which.plots){
      plots$confusion <- heatmapFromMatrix(
         tab, show.labels=T, label.size=label.size,
         x.lab='Actual class', y.lab='Predicted class', invert.y=T,
         legend.name=if(rel.values){ 'Column fraction' } else { 'Counts' }
      ) +
         theme(
            axis.text.x.bottom=element_blank(),
            #axis.ticks=element_line(),
            legend.position='none'
         )
   }

   ## Correct/incorrect samples per class ----------------------------
   ##
   class_counts <- table(actual)
   class_counts <- class_counts[class_order]
   class_counts <- c(Total=length(actual), class_counts)

   correct_counts <- unlist(lapply(split(data.frame(actual, predicted), actual), function(i){
      sum(i$actual==i$predicted)
   }))
   correct_counts <- correct_counts[class_order]
   correct_counts <- c(Total=sum(actual==predicted),correct_counts)

   counts <- data.frame(Total=class_counts, Correct=correct_counts)
   counts$Incorrect <- counts$Total - counts$Correct

   ##
   if('counts' %in% which.plots){
      pd_counts <- structure(
         reshape2::melt(as.matrix(counts)),
         names=c('class','measure','value')
      )

      pd_counts$class <- as.character(pd_counts$class)
      pd_counts$Total <- counts[pd_counts$class,'Total']
      pd_counts$frac <- pd_counts$value / pd_counts$Total
      pd_counts$frac[pd_counts$measure=='Total'] <- NA
      pd_counts$frac[pd_counts$measure=='Incorrect'] <- -pd_counts$frac[pd_counts$measure=='Incorrect']

      pd_counts$class[pd_counts$class=='Total'] <- 'OVERALL'
      pd_counts$class <- factor(pd_counts$class, unique(pd_counts$class))

      plots$counts <- ggplot(pd_counts, aes(y=measure, x=class)) +
         geom_tile(aes(fill=frac), color='grey') +
         geom_text(aes(label=value), size=label.size) +
         geom_hline(yintercept=1.5, size=0.5) +
         scale_fill_distiller(palette='RdYlGn', na.value='white', limits=c(-1,1), direction=1, guide=F) +
         scale_x_discrete(expand=c(0,0), position='top') +
         scale_y_discrete(expand=c(0,0)) +
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            axis.text.x.bottom=element_text(angle=90, vjust=0.5, hjust=1),
            axis.text.x.top=element_text(angle=90, vjust=0.5, hjust=0),
            axis.title=element_blank(),
            axis.ticks=element_blank()
         )
   }

   ## Performance stats ----------------------------
   if('perf' %in% which.plots){
      confusion <- mltoolkit::confusionMatrix(
         predicted=mltoolkit::oneHotEncode(predicted),
         actual=actual,
         simplify=T
      )

      #metrics=c('f1','precision')
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

      ## Capitalize first letter of metric names
      firstupper <- function(x) {
         substr(x, 1, 1) <- toupper(substr(x, 1, 1))
         x
      }
      colnames(perf) <- firstupper(colnames(perf))

      ## Plot
      plots$perf <- heatmapFromMatrix(
         t(perf),
         palette='RdYlGn', palette.direction=1,
         show.labels=T, label.size=label.size,
         invert.y=T, y.lab='Perf.', legend.name='Metric value', x.title.position='top'
      ) +
         guides(fill=F) +
         theme(
            axis.text.x.bottom=element_blank(),
            axis.title=element_blank()
         )
   }

   ## Combine --------------------------------
   plots <- plots[which.plots]

   ## Assign rel heights based on `rel.heights` names or index
   if(length(names(rel.heights))!=0){
      rel_heights <- rel.heights[names(plots)]
   } else {
      rel_heights <- rel.heights[1:length(plots)]
      names(rel_heights) <- names(plots)
   }

   ##
   plots <- lapply(plots, function(i){
      i + theme(
         axis.text.x.bottom=element_blank(),
         axis.text.x.top=element_blank(),
         axis.ticks.y=element_line()
      )
   })

   ## Add x-axis labels and/or title to top plot
   plots[[1]] <- plots[[1]] +
      theme(
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.ticks.x=element_line()
      )
   if(!is.null(plot.title)){
      plots[[1]] <- plots[[1]] + ggtitle(plot.title)
   }

   ## Add x-axis labels bottom plot
   plots[[length(plots)]] <- plots[[length(plots)]] +
      theme(
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x.bottom=element_text(),
         axis.ticks.x=element_line()
      )

   cowplot::plot_grid(
      plotlist=plots,
      ncol=1, align='v', axis='tblr', rel_heights=rel_heights
   )

}

####################################################################################################
#' Calculate the fraction of correctly predicted samples
#'
#' @param actual A factor of the actual classes
#' @param predicted A factor of the predicted classes
#' @param rm.non.pred.classes Remove stats for classes that are not in the factor levels of
#' `predicted`
#'
#' @return A dataframe
#' @export
#'
calcFracCorrect <- function(actual, predicted, rm.non.pred.classes=T){

   if(!is.factor(actual)){ stop('`actual` must be a factor') }
   if(!is.factor(predicted)){ stop('`predicted` must be a factor') }

   ## Format factor levels
   pred_levels <- levels(predicted)

   uniq_responses <- sort(unique(
      c(
         as.character(actual),
         as.character(predicted)
      )
   ))

   actual <- factor(actual, uniq_responses)
   predicted <- factor(predicted, uniq_responses)

   n_correct <- diag(table(predicted, actual))
   n_total <- table(actual)

   ## Main
   out <- data.frame(
      #group=i,
      class=names(n_correct),
      n_correct=as.integer(n_correct),
      n_total=as.integer(n_total)
   )

   if(rm.non.pred.classes){
      out <- out[out$class %in% pred_levels,]
   }

   out <- rbind(
      data.frame(
         class='All',
         n_correct=sum(out$n_correct),
         n_total=sum(out$n_total)
      ),
      out
   )

   out$frac_correct <- out$n_correct / out$n_total

   return(out)
}




