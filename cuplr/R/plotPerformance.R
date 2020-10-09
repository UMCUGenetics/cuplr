####################################################################################################
#' Apply a summary function to a list of matrixes of the same dimensions
#'
#' @param l A list of matrices
#' @param func A summary function (default is mean())
#' @param as.matrix If TRUE, will return a single matrix with the aggregated values. If FALSE, will
#' return a melted dataframe.
#' @param value.name A custom column name for the value column of the melted dataframe
#'
#' @export
#'
aggregateMatrixList <- function(l, func=function(x){ mean(x, na.rm=T) }, as.matrix=F, value.name=NULL){
   #l=lapply(cv_out,`[[`,'imp')

   all_cols <- unique(unlist(lapply(l, colnames)))
   l_melt <- lapply(l, function(i){
      #i=l[[3]]
      df <- as.data.frame(i)
      missing_cols <- all_cols[!(all_cols %in% colnames(df))]
      df[,missing_cols] <- NA
      df_melt <- reshape2::melt(as.matrix(df))
      colnames(df_melt) <- c('class','feature','value')

      return(df_melt)
   })

   df_pre <- suppressWarnings({
      Reduce(function(x,y){
         merge(x,y,by=c(1,2), all=T, sort=F)
      }, l_melt)
   })

   df <- df_pre[,c(1,2)]
   m_values <- as.matrix(df_pre[,-c(1,2)])
   df$value <- apply(m_values,1,function(i){ func(i) })



   #m_values <- do.call(cbind, lapply(l_melt,`[[`,'value'))
   #df <- l_melt[[1]]
   #df$value <- apply(m_values,1,function(i){ func(i) })

   if(as.matrix){
      out <- reshape2::acast(df, class ~ feature)
      out <- out[,colnames(l[[1]])]
      return(out)
   } else {
      if(!is.null(value.name)){ colnames(df)[3] <- value.name }
      return(df)
   }
}

####################################################################################################
#' Make barplots from a matrix of feature importances (for multiclass classification)
#'
#' @param m A numeric matrix where columns represent the features and rows represent the labels
#' @param cv_out Alternative input to `m`. Cross validation output in the form of a list that has
#' 'imp' object (e.g. cv_out[[1]]$imp)
#' @param top.n Top number of features to show
#' @param infer.feature.type Determine the feature type based on the tag/prefix. Everything before
#' the first dot (.) is considered the feature tag
#' @param n.row Number of facet rows
#' @param n.col Number of facet columns
#' @param feature.type.colors A character vector of color hex codes with names being the feature
#' types
#'
#' @return A ggplot object
#' @export
#'
plotTopFeatures <- function(
   m=NULL, cv_out=NULL, top.n=10, infer.feature.type=F,
   n.row=NULL, n.col=NULL,
   feature.type.colors=NULL
){
   if(!is.null(cv_out)){
      m <- aggregateMatrixList(lapply(cv_out,`[[`,'imp'), as.matrix=T)
   }

   df <- do.call(rbind, lapply(rownames(m), function(i){
      v <- sort(m[i,], decreasing=T)[1:top.n]
      data.frame(class=i, feature=names(v), value=v, index=1:top.n, row.names=NULL)
   }))
   df <- forceDfOrder(df)

   if(infer.feature.type){
      df$feature_type <- factor(
         gsub('[.].+$','',df$feature),
         levels=unique(gsub('[.].+$','',colnames(m)))
      )
   } else {
      df$feature_type <- 'none'
   }

   df$label <- as.character(df$feature)
   df$label[df$value<=0] <- ''

   label_y_pos <- max(df$value) * 0.05

   if(is.null(feature.type.colors)){
      color_pal <- c(
         RColorBrewer::brewer.pal(12, 'Set3'),
         RColorBrewer::brewer.pal(9, 'Pastel1'),
         RColorBrewer::brewer.pal(8, 'Pastel2')
      )
      color_pal <- structure(
         color_pal[1:length(levels(df$feature_type))],
         names=levels(df$feature_type)
      )
   } else if(feature.type.colors=='auto'){
      color_pal <- c(
         sigs="#8DD3C7",
         kataegis="#FFFFB3",
         gene_amp="#FB8072",
         gene_def="#80B1D3",
         aneuploidy="#FDB462",
         #arm_loss="#FDB462",
         #arm_gain="#B3DE69",
         purple='#BC80BD',
         sv_types="#FCCDE5",
         viral_ins="#D9D9D9",
         fusion="#FFED6F",
         rep_elem="#CCEBC5",
         rmd="#BEBADA"
      )
   } else {
      color_pal <- feature.type.colors
   }


   require(ggplot2)
   p <- ggplot(df, aes(x=index, y=value)) +
      facet_wrap(~class, nrow=n.row, ncol=n.col) +

      #geom_bar(stat='identity', fill='#659F5B') +
      geom_bar(aes(fill=feature_type), stat='identity') +
      #scale_fill_brewer(palette='Set3') +
      scale_fill_manual(values=color_pal, limits=names(color_pal)) +

      geom_text(aes(label=label, y=label_y_pos), angle=90, hjust=0, size=2.5) +

      labs(y='Feature importance', x='Index', fill='Feature type') +

      theme(
         panel.grid.minor=element_blank()
      )

   if(length(unique(df$feature_type))==1){
      p <- p + guides(fill=F)
   }

   return(p)
}

####################################################################################################
#' Make heatmap from a matrix of feature importances (for multiclass classification)
#'
#' @param m A numeric matrix where columns represent the features and rows represent the labels
#' @param cv_out Alternative input to `m`. Cross validation output in the form of a list that has
#' 'imp' object (e.g. cv_out[[1]]$imp)
#' @param top.n Top number of features to show
#' @param min.imp Features below this importance value are excluded from the plot
#' @param invert.y Invert the y-axis?
#'  @param sort.features Order features alphabetically?
#'
#' @return A ggplot object
#' @export
#'
plotFeatureImpHeatmap <- function(
   m=NULL, cv_out=NULL, top.n=50, min.imp=NULL, invert.y=T, sort.features=F
){
   if(!is.null(cv_out)){
      m <- aggregateMatrixList(lapply(cv_out,`[[`,'imp'), as.matrix=T)
   }
   if(sort.features){
      m <- m[,order(colnames(m))]
   }

   df <- reshape2::melt(as.matrix(m))
   colnames(df) <- c('class','feature','value')

   if(!is.null(top.n)){
      df2 <- df[order(df$value, decreasing=T),]
      feature_whitelist <- unique(df2$feature)[1:top.n]
      df <- df[df$feature %in% feature_whitelist,]
   } else if(!is.null(min.imp)){
      feature_whitelist <- as.character(unique(df[df$value >= min.imp,'feature']))
      df <- df[df$feature %in% feature_whitelist,]
   }

   df <- forceDfOrder(df)
   df$index <- as.integer(df$feature)

   if(invert.y){
      df$class <- factor(df$class, rev(levels(df$class)))
   }

   n_sel_features <- length(unique(df$feature))
   xlabel <-
      if(ncol(m) < n_sel_features ){
         'Features'
      } else {
         paste0('Features (top ',n_sel_features,'/',ncol(m),')')
      }

   require(ggplot2)
   p <- ggplot(df, aes(y=class,x=index)) +
      geom_tile(aes(fill=value)) +
      scale_x_continuous(
         sec.axis=dup_axis(), breaks=unique(df$index), labels=levels(df$feature),
         expand=c(0,0)
      ) +
      scale_fill_distiller(palette='YlGnBu') +

      labs(y='Class',x=xlabel,fill='Feat.\nimp.') +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5)
      )

   if(!is.null(min.imp)){
      p <- p + xlab(sprintf('Features\n(min.imp>=%s in at least 1 class)', min.imp))
   }

   return(p)
}

####################################################################################################
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
   rel.heights=c(0.3, 0.05, 1),

   ## Confusion heatmap
   sort=F, rel.values=T,

   ## Performance metrics
   metrics=c('f1','prec','tpr'), show.weighted.mean=T
){
   require(ggplot2)
   if(!is.null(cv_out)){
      actual <- unlist(lapply(cv_out,function(i){ i$test_set$actual }))
      predicted <- unlist(lapply(cv_out,function(i){ i$test_set$predicted }))
   }

   ## Large confusion matrix ----------------------------------------------------------------
   tab <- table(predicted, actual)

   if(rel.values){
      tab <- apply(tab,2,function(i){ i/sum(i) })
      tab <- round(tab,2)
   }

   class_order <- if(sort){ names(sort(diag(tab), decreasing=T)) } else { names(diag(tab)) }
   tab <- tab[class_order,class_order]

   tab <- cbind(NA,tab) ## Add dummy column for mean performance

   p_heatmap <- plotHeatmapFromMatrix(
      tab, show.labels=T, x.lab='Actual\n(Columns: prop. misclassified as which class)', y.lab='Predicted', invert.y=T,
      #palette=if(rel.values){ 'YlGnBu' } else { 'none' },
      legend.name=if(rel.values){ 'Column fraction' } else { 'Counts' }
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

   p_counts <-
      plotHeatmapFromMatrix(t(as.matrix(counts)), palette='none', show.labels=T) +
      geom_hline(yintercept=1.5, size=0.5) +
      theme(
         #panel.grid=element_blank(),
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

   perf <- perf[class_order,]

   ## Summary stats
   if(show.weighted.mean){
      perf_summary <- rbind(
         OVERALL=apply(perf,2,function(i){
            weights <- table(actual)/length(actual)
            weights <- weights[names(i)]
            weighted.mean(i, weights)
         })
      )
   } else {
      perf_summary <- rbind(
         OVERALL = apply(perf,2,mean)
      )
   }
   perf <- rbind(perf_summary, perf)
   perf <- round(perf, 2)

   p_perf <- plotHeatmapFromMatrix(
         t(perf), palette='RdYlGn', palette.direction=1, show.labels=T,
         y.lab='Perf.', legend.name='Metric value', x.title.position='top'
      ) +
      theme(
         axis.title.x=element_blank()
      )

   ## Combine --------------------------------
   cowplot::plot_grid(
      p_perf, p_counts, p_heatmap,
      ncol=1, align='v', axis='tblr', rel_heights=c(0.3, 0.15, 1)
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

####################################################################################################
mkTrainingReport <- function(dir, verbose=T){
   if(F){
      dir='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.08a_blacklistSamples_probWeighRf/'
      verbose=T
   }

   plots_dir <- paste0(dir,'/plots/')
   dir.create(plots_dir, showWarnings=F)

   if(verbose){ message('Loading feat imp data...') }
   imp <- readRDS(paste0(dir,'/imp.rds'))
   m_imp <- aggregateMatrixList(imp, as.matrix=T)

   if(verbose){ message('Plotting imp barplots...') }
   pdf(paste0(plots_dir,'/imp_barplots.pdf'), 16, 10)
   plot( plotTopFeatures(m_imp, top.n=40, infer.feature.type=T, n.col=4) )
   dev.off()

   if(verbose){ message('Loading feat imp heatmap...') }
   pdf(paste0(plots_dir,'/imp_heatmap.pdf'), 17, 8)
   plot( plotFeatureImpHeatmap(m_imp, top.n=150) )
   dev.off()
   #plotFeatureImpHeatmap(m_imp, min.imp=1)

   if(verbose){ message('Loading test set data...') }
   test_set <- readRDS(paste0(dir,'/test_set.rds'))
   actual <- unlist(lapply(test_set,`[[`,'actual'))
   predicted <- unlist(lapply(test_set,`[[`,'predicted'))
   prob <- as.data.frame(do.call(rbind, lapply(test_set, function(i){ i$probabilities })))

   if(verbose){ message('Plotting perf heatmap...') }
   pdf(paste0(plots_dir,'/perf_heatmap.pdf'), 11, 8.5)
   suppressWarnings({
      plot( plotPerfHeatmap(actual, predicted, show.weighted.mean=T) )
   })
   dev.off()

   if(verbose){ message('Plotting false negative rates...') }
   pdf(paste0(plots_dir,'/fnr_heatmap.pdf'), 11, 8)
   plot( plotFnr(actual, prob) )
   dev.off()
}

