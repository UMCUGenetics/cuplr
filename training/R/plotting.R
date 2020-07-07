####################################################################################################
aggregateMatrixList <- function(l, func=function(x){ mean(x, na.rm=T) }, value.name=NULL, as.matrix=F){
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
      return(reshape2::acast(df, class ~ feature))
   } else {
      if(!is.null(value.name)){ colnames(df)[3] <- value.name }
      return(df)
   }
}

####################################################################################################
plotTopFeatures <- function(
   m=NULL, cv_out=NULL, top.n=10, infer.feature.type=F,
   nrow=NULL, ncol=NULL
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
   
   label_y_pos <- max(df$value) * 0.05
   
   require(ggplot2)
   p <- ggplot(df, aes(x=index, y=value)) +
      facet_wrap(~class, nrow=nrow, ncol=ncol) +
      
      #geom_bar(stat='identity', fill='#659F5B') +
      geom_bar(aes(fill=feature_type), stat='identity') +
      scale_fill_brewer(palette='Set3') +
      
      geom_text(aes(label=feature, y=label_y_pos), angle=90, hjust=0, size=2.5) +
      
      labs(y='Feature importance', x='Index', fill='Feature type') +
      
      theme(
         panel.grid.minor=element_blank()
      )
   
   if(length(unique(df$feature_type))==1){
      p <- p + guides(fill=F)
   }
   
   return(p)
}

plotFeatureImpHeatmap <- function(
   m=NULL, cv_out=NULL, min.imp=NULL
){
   if(!is.null(cv_out)){
      m <- aggregateMatrixList(lapply(cv_out,`[[`,'imp'), as.matrix=T)
   }
   
   df <- reshape2::melt(as.matrix(m))
   colnames(df) <- c('class','feature','value')
   
   if(!is.null(min.imp)){
      feature_whitelist <- as.character(unique(df[df$value >= min.imp,'feature']))
      df <- df[df$feature %in% feature_whitelist,]
   }
   
   df <- forceDfOrder(df)
   
   require(ggplot2)
   p <- ggplot(df, aes(y=class,x=feature)) + 
      geom_tile(aes(fill=value)) +
      scale_fill_distiller(palette='YlGnBu') +
      
      labs(y='Class',x='Features',fill='Imp.') +
      
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
      )
   
   if(!is.null(min.imp)){
      p <- p + xlab(sprintf('Features\n(min.imp>=%s in at least 1 class)', min.imp))
   }
   
   return(p)
}

####################################################################################################
plotHeatmapFromMatrix <- function(
   m, x.lab=NULL, y.lab=NULL, legend.name='value', 
   show.labels=F, palette='YlGnBu', palette.direction=-1, invert.y=F
){
   
   #m=t(imp)
   if(!is.matrix(m) & !is.numeric(m)){ stop('`m` must be a numeric matrix') }
   m_melt <- reshape2::melt(m)
   colnames(m_melt)[1:2] <- c('y','x')
   
   if(invert.y){
      m_melt$y <- factor(m_melt$y, rev(levels(m_melt$y)))
   }
   
   p <- ggplot(m_melt, aes(x, y, fill=value)) +
      geom_tile(color='grey') +
      scale_fill_distiller(
         name=legend.name, 
         palette=palette, direction=palette.direction, 
         guide=guide_colorbar(
            frame.colour='black', ticks.colour='black', 
            direction='horizontal', title.position='top', reverse=T, barwidth=4, barheight=1,
            label.theme=element_text(angle=90, vjust=0.5, size=10)
         )
      ) +
      #scale_x_discrete(expand=c(0,0)) +
      #scale_y_discrete(expand=c(0,0)) + 
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x.bottom=element_text(angle=90, vjust=0.5, hjust=1),
         axis.text.x.top=element_text(angle=90, vjust=0.5, hjust=0),
         axis.ticks=element_blank()
      )
   
   if(show.labels){ p <- p + geom_text(aes(label=value), size=3.5) }
   if(!is.null(x.lab)){ p <- p + xlab(x.lab) }
   if(!is.null(y.lab)){ p <- p + ylab(y.lab) }
   
   return(p)
}

plotPerfHeatmap <- function(
   actual=NULL, predicted=NULL, cv_out=NULL, 
   rel.heights=c(0.3, 1),
   
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
   tab <- table(actual, predicted)
   
   if(rel.values){
      tab <- t(apply(tab,1,function(i){ i/sum(i) }))
      tab <- round(tab,2)
   }
   
   class_order <- if(sort){ names(sort(diag(tab), decreasing=T)) } else { names(diag(tab)) }
   tab <- tab[class_order,class_order]
   
   tab <- cbind(NA,tab) ## Add dummy column for mean performance
   
   p_heatmap <- plotHeatmapFromMatrix(
      tab, show.labels=T, y.lab='Actual', x.lab='Predicted', invert.y=T,
      legend.name=if(rel.values){ 'Prop. predicted' } else { 'Counts' }
   )
   
   ## Samples per class ----------------------------------------------------------------
   class_counts <- table(actual)
   class_counts <- class_counts[colnames(tab)]
   
   class_counts[1] <- length(actual)
   names(class_counts)[1] <- 'total'
   
   class_counts <- data.frame(class=names(class_counts), counts=as.numeric(class_counts))
   class_counts <- forceDfOrder(class_counts)
   
   p_class_counts <- ggplot(class_counts, aes(x=class, y=1, label=counts)) +
      geom_tile(fill='white',color='black') +
      geom_text(size=3.5) +
      ylab('# samples') +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_text(angle=0, vjust=0.5, size=10),
         axis.ticks=element_blank()
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
         WEIGHTED_MEAN=apply(perf,2,function(i){ 
            weights <- table(actual)/length(actual)
            weights <- weights[names(i)]
            weighted.mean(i, weights) 
         })
      )
   } else {
      perf_summary <- rbind(
         MEAN = apply(perf,2,mean)
      )
   }
   perf <- rbind(perf_summary, perf)
   perf <- round(perf, 2)
   
   p_perf <- plotHeatmapFromMatrix(
         t(perf), palette='RdYlGn', palette.direction=1, show.labels=T, 
         y.lab='Perf.', legend.name='Metric value'
      ) +
      scale_x_discrete(position='top') +
      theme(
         axis.title.x=element_blank()
      )
   
   ## Combine --------------------------------
   cowplot::plot_grid(
      p_perf, p_class_counts, p_heatmap, 
      ncol=1, align='v', axis='tblr', rel_heights=c(0.3, 0.05, 1)
   )
   
}

####################################################################################################
# plotPerfCustom <- function(
#    actual=NULL, predicted=NULL, cv_out=NULL, metrics=c('f1','prec','rec'), 
#    sort=T, show.summary.stats=T
# ){
#    #metrics=c('f1','prec','rec')
#    
#    if(!is.null(cv_out)){
#       actual <- unlist(lapply(cv_out,function(i){ i$test_set$actual }))
#       predicted <- unlist(lapply(cv_out,function(i){ i$test_set$predicted }))
#    }
#    
#    confusion <- mltoolkit::confusionMatrix(
#       predicted=mltoolkit::oneHotEncode(predicted), 
#       actual=actual, 
#       simplify=T
#    )
#    
#    perf <- mltoolkit::calcPerf(confusion, metrics)
#    rownames(perf) <- perf[,1]; perf[,1] <- NULL
#    
#    if(sort){
#       perf <- perf[order(perf[,1], decreasing=T),]
#    }
#    
#    if(show.summary.stats){
#       perf_agg <- rbind(
#          MEAN = apply(perf,2,mean),
#          WEIGHTED_MEAN = apply(perf,2,function(i){ 
#             weights <- table(actual)/length(actual)
#             weights <- weights[names(i)]
#             weighted.mean(i, weights) 
#          })
#       )
#       
#       perf <- rbind(perf_agg,perf)
#    }
#    
#    perf <- cbind(class=rownames(perf), perf)
#    rownames(perf) <- NULL
#    
#    pd <- reshape2::melt(perf,'class')
#    colnames(pd) <- c('class','metric','value')
#    pd <- forceDfOrder(pd)
#    
#    ggplot(pd, aes(x=class, y=value)) +
#       facet_wrap(~metric, ncol=1) +
#       geom_bar(stat='identity', fill='#d26c6c') +
#       geom_text(aes(label=round(value,2), y=0.02), angle=90, hjust=0, vjust=0.5) +
#       theme_bw() +
#       theme(
#          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
#       )
#    
# }










