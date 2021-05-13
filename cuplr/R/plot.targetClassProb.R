#' Plot probability of target class
#'
#' @param actual A vector of the actual classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param report A list with the objects with the names: prob, class_actual
#' @param output If 'plot', returns a ggplot object, else 'raw' returns the raw data
#' @param cutoffs A numeric vector of cutoffs. Applies only when `plot.type`='per_cutoff'
#' @param force.discrete.x.axis Force discrete x axis
#' @param show.labels Show labels indicating the numerical values at each point?
#' @param show.bars For each sample, show bars colored by the predicted class? Only applies when
#' `output`=='plot.sorted_probs'
#' @param show.all.classes.stats Show the panel with the stats for all samples (i.e. not split by
#' actual class)
#'
#' @description For samples within each actual class, plot the probability of said actual class.
#' The type of plot returned depends on `output`:
#'   'plot.sorted_probs'  x: cutoff, y: probability of correct prediction
#'   'plot.per_cutoff'    x: cutoff, y: fraction of samples >= cutoff
#'
#' @return A ggplot object
#' @export
#'
targetClassProb <- function(
   actual=NULL, probs=NULL, report=NULL,
   output=c('plot.sorted_probs','plot.per_cutoff','raw.sorted_probs','raw.per_cutoff'),
   show.labels=F, show.bars=F,
   cutoffs=seq(0, 1, 0.1), force.discrete.x.axis=F,
   show.all.classes.stats=if(show.bars){ F } else { T }
){
   if(F){
      report <- pred_reports$holdout
      report <- pred_reports$CUP

      actual=report$class_actual
      #actual=report$class_pred
      #predicted=report$class_pred
      probs=report$prob_scaled
      force.discrete.x.axis=F
      output=c('plot.sorted_probs','plot.per_cutoff','raw.sorted_probs','raw.per_cutoff')
      cutoffs=seq(0, 1, 0.1)
      show.all.classes.stats=T
   }

   ## Init --------------------------------
   require(ggplot2)

   if(!is.null(report)){
      actual <- report$class_actual
      probs <- report$prob
   }
   predicted <- factor(
      colnames(probs)[ max.col(probs) ],
      colnames(probs)
   )

   if(length(actual)!=nrow(probs)){
      stop('length(actual) and nrow(probs) must be equal')
   }

   if(!is.factor(actual)){
      warning('`actual` is not a factor. Factor levels will be set automatically')
      actual <- as.factor(actual)
   }

   output <- match.arg(output, c('plot.sorted_probs','plot.per_cutoff','raw.sorted_probs','raw.per_cutoff'))

   ## Get prob of target cancer type --------------------------------
   df <- data.frame(actual,predicted,row.names=NULL)
   df <- cbind(sample=rownames(probs), df)

   probs_melt <- reshape2::melt(probs)
   colnames(probs_melt) <- c('sample','binary_rf','prob')

   df$prob_target_class <- probs_melt$prob[
      match(
         paste0(df$sample, as.character(df$actual)),
         paste0(probs_melt$sample, probs_melt$binary_rf)
      )
   ]
   rm(probs_melt)

   uniq_classes <- levels(actual)

   if(show.all.classes.stats){
      df <- rbind(
         within(df,{ actual <- 'All' }),
         df
      )
      df$actual <- factor(df$actual, c('All',uniq_classes))
   }

   ## Sorted probs --------------------------------
   df <- df[order(df$actual, -df$prob_target_class),]
   df <- do.call(rbind, lapply(split(df, df$actual), function(i){
      #i=split(df, df$actual)[[1]]
      if(nrow(i)>0){
         i$index <- nrow(i):1
      } else {
         i$index <- integer()
      }
      return(i)
   })); rownames(df) <- NULL

   n_samples <- table(df$actual)
   df$n_samples <- as.integer( n_samples[ df$actual ] )
   df$actual.n_samples <- with(df,paste0(
      actual,' (',n_samples,')'
   ))

   if(output=='raw.sorted_probs'){ return(df) }

   if(output=='plot.sorted_probs'){

      p <- ggplot(df, aes(x=index, y=prob_target_class, group=1)) +
         facet_wrap(~actual.n_samples, scales='free_x') +
         scale_y_continuous(name='Prob. of target class', limits=c(0,1)) +
         scale_x_continuous(name='Sample rank') +
         theme_bw() +
         theme(
            panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank()
         )

      if(show.bars){
         color_pal <- c(
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", ## Set3
            "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", ## Set2
            "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", ## Pastel1
            "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999" ## Set1
         )
         class_colors <- structure( color_pal[1:length(uniq_classes)], names=uniq_classes )

         p <- p +
            geom_bar(aes(fill=predicted), stat='identity', width=1) +
            scale_fill_manual(name='Predicted class', values=class_colors)
      }

      ## Outline
      p <- p +
         geom_step(aes(x=index-0.5), size=0.3) +
         geom_segment(
            data=df[df$index==df$n_samples,],
            aes(x=index-0.5, xend=index+0.5, y=prob_target_class, yend=prob_target_class),
            size=0.3
         )

      return(p)
   }

   ## Per cutoff --------------------------------
   ## No. samples above cutoff
   m_agg <- do.call(cbind, lapply(cutoffs, function(i){
      #i=cutoffs[1]
      actual_classes <- df$actual[df$prob_target_class>=i]
      table(actual_classes)
   }))
   colnames(m_agg) <- cutoffs
   df_agg <- reshape2::melt(m_agg)
   colnames(df_agg) <- c('actual','cutoff','n_ge_cutoff')
   df_agg <- df_agg[order(df_agg$actual),]

   ## Prop above cutoff
   df_agg$n_samples <- as.integer( n_samples[ df_agg$actual ] )
   df_agg$frac_ge_cutoff <- df_agg$n_ge_cutoff / df_agg$n_samples
   df_agg <- df_agg[!is.na(df_agg$frac_ge_cutoff),]

   if(output=='raw.per_cutoff'){ return(df_agg) }

   ## Labels
   df_agg$actual.n_samples <- with(df_agg,paste0(
      actual,' (',n_samples,')'
   ))

   df_agg$label <- with(df_agg, paste0(
      round(frac_ge_cutoff,2),' (',n_ge_cutoff,')'
   ))
   df_agg$label.hjust <- 1
   df_agg$label.hjust[df_agg$frac_ge_cutoff<0.5] <- 0
   df_agg$label.ypos <- df_agg$frac_ge_cutoff + ifelse(df_agg$frac_ge_cutoff<0.5, 0.05, -0.05)

   if(force.discrete.x.axis){
      df_agg$cutoff <- factor(df$cutoff, cutoffs)
   }

   p <- ggplot(df_agg, aes(x=cutoff, y=frac_ge_cutoff)) +
      facet_wrap(~actual.n_samples) +
      geom_point(color='#6F9DC5') +
      geom_line(aes(group=1), color='#6F9DC5') +
      labs(y='Frac. samples >= cutoff', x='Cutoff') +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         panel.grid.minor.y=element_blank()
      )

   if(show.labels){
      p <- p + geom_text(aes(label=label, y=label.ypos, hjust=label.hjust), angle=90, size=2.5)
   }

   return(p)
}



