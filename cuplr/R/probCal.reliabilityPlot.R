#' Reliability plot
#'
#' @description Bin probabilities and calculate the proportion of positive class samples per bin.
#' The proportion of positive class samples (y-values) will be plotted against the midpoint of each
#' probability bin (x-values), and a curve will connect the points. A 'filtered' curve is also shown
#' where points supported by too few samples (<`min.samples.per.bin`) will be ignored for plotting
#' this curve.
#'
#' @param actual A vector of the actual classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param report A list with the objects with the names: prob, class_actual
#' @param n.bins Number of (equidistant bins)
#' @param bootstrap.iters Number of bootstrap iterations to calculate error bars (standard deviation)
#' @param bootstrap.prop Proportion of all samples to sample from for boostrapping
#' @param seed Bootstrap seed
#' @param show.labels Show labels at each point
#' @param min.samples.per.bin Points with fewer samples than this will no be considered for the
#' overlaid (filtered) reliability curve
#' @param verbose Show progress messages?
#'
#' @return A ggplot object
#' @export
#'
reliabilityPlot <- function(
   actual=NULL, probs=NULL, report=NULL,
   n.bins=10, bootstrap.iters=50, bootstrap.prop=0.66, seed=1,
   show.labels=T, min.samples.per.bin=4, verbose=F
){

   # if(F){
   #    report=pred_reports$CV
   #
   #    actual=report$class_actual
   #    probs=report$prob_scaled
   #    #probs=report$prob
   #
   #    n.bins=10
   #    bootstrap.iters=50
   #    bootstrap.prop=0.66
   #    seed=1
   #    show.labels=T
   #    min.samples.per.bin=3
   # }

   ## Init --------------------------------
   if(!is.null(report)){
      actual <- report$class_actual
      probs <- report$prob
   }

   if(length(actual)!=nrow(probs)){
      stop('length(actual) and nrow(probs) must be equal')
   }

   if(!all(colnames(probs) %in% unique(actual))){
      stop('`colnames(probs)` must have the same classes as in `actual`')
   }

   ## Calc fraction from positive class --------------------------------
   ## Bins
   bins <- seq(from=0, to=1, by=1/n.bins)
   bin_mids <- (
      bins[2:length(bins)]
      + bins[-length(bins)]
   ) / 2

   ## Bootstrap init
   n_samples <- length(actual)
   bootstrap_size <- round( bootstrap.prop * n_samples )

   set.seed(seed)
   uniq_classes <- colnames(probs)
   if(verbose){ message('Calculating prop. of positive class samples per bin...') }
   stats <- lapply(uniq_classes, function(i){
      #i='Lung_SC'

      if(verbose){ message('> ',i) }

      ## Prep data
      df <- data.frame(
         prob=probs[,i],
         response=actual==i,
         row.names=NULL
      )

      df$bin_mid <- cut(df$prob, breaks=bins, include.lowest=T, labels=bin_mids)

      ## Main
      calcStats <- function(df){
         #df_split <- split(bootstrap_out, bootstrap_out$bin_mid)
         df_split <- split(df, df$bin_mid)
         agg <- do.call(rbind, lapply(df_split, function(j){
            #j=df_split[[1]]
            if(nrow(j)!=0){
               data.frame(
                  prob.mean=mean(j$prob, na.rm=T),
                  frac_pos=mean(j$response, na.rm=T),
                  #frac_pos.sd=sd(j$response, na.rm=T),
                  n_pos=sum(j$response, na.rm=T),
                  n_samples=nrow(j)
               )
            } else {
               data.frame(
                  prob.mean=NA,
                  frac_pos=NA,
                  #frac_pos.sd=NA,
                  n_pos=0,
                  n_samples=0
               )
            }

         }))
         #agg$frac_pos.sd[ is.na(agg$frac_pos.sd) ] <- 0
         agg <- cbind(bin_mid=as.numeric(rownames(agg)), agg)
         return(agg)
      }

      # if(is.null(bootstrap.iters)){
      #    return( calcStats(df) )
      # }

      bootstrap_out <- do.call(rbind, lapply(1:bootstrap.iters, function(j){
         res <- calcStats( df[sample(1:n_samples, bootstrap_size),] )
         res$boot_iter <- j
         return(res)
      }))
      rownames(bootstrap_out) <- NULL

      bootstrap_out <- split(bootstrap_out, bootstrap_out$bin_mid)

      out <- do.call(rbind, lapply(bootstrap_out, function(j){
         #j=bootstrap_out[[1]]

         if(nrow(j)!=0){
            data.frame(
               prob.mean=mean(j$prob.mean, na.rm=T),
               frac_pos=mean(j$frac_pos, na.rm=T),
               frac_pos.sd=sd(j$frac_pos, na.rm=T),
               n_pos.mean=mean(j$n_pos, na.rm=T),
               n_samples.mean=mean(j$n_samples, na.rm=T)
            )
         } else {
            data.frame(
               prob.mean=NA,
               frac_pos=NA,
               frac_pos.sd=NA,
               n_pos.mean=0,
               n_samples.mean=0
            )
         }
      }))
      out$frac_pos.sd
      out <- data.frame(bin_mid=as.numeric(rownames(out)), out, row.names=NULL)

      return(out)
   })
   names(stats) <- uniq_classes

   # ## Make curve isotonic (y always increasing) --------------------------------
   # stats <- lapply(stats, function(i){
   #    #i=stats[[1]]
   #
   #    current_frac_pos <- i$frac_pos[1]
   #    #current_frac_pos=0.892015873
   #    i$frac_pos_iso <- sapply(i$frac_pos, function(j){
   #       #j=0.608247423
   #       if(j > current_frac_pos){
   #          current_frac_pos <<- j
   #       }
   #       return(current_frac_pos)
   #    })
   #
   #    return(i)
   # })

   stats <- do.call(rbind, lapply(names(stats), function(i){
      cbind(class=i, stats[[i]])
   }))
   stats$class <- factor(stats$class, uniq_classes)

   ## Plot --------------------------------
   if(verbose){ message('Generating reliability plot...') }
   pd <- stats
   pd <- pd[!is.na(pd$frac_pos),]

   pd$n_pos <- round(pd$n_pos.mean)
   pd$n_samples <- round(pd$n_samples.mean)
   pd$ge_min_samples_per_bin <- pd$n_samples >= min.samples.per.bin

   pd$label <- paste0(pd$n_pos,'/',pd$n_samples)
   pd$label.ypos <- 0
   pd$label.ypos[pd$bin_mid<0.5] <- 1

   pd$label.hjust <- 0
   pd$label.hjust[pd$bin_mid<0.5] <- 1

   color_raw <- 'grey50'
   color_filt <- 'indianred'

   p <- ggplot(pd, aes(x=bin_mid, y=frac_pos)) +
      facet_wrap(class~.) +

      geom_abline(slope=1, intercept=0, linetype='dashed', color='lightgrey') +

      ## Original curve
      geom_path(size=0.3, color=color_raw) +
      geom_linerange(
         aes(ymin=frac_pos-frac_pos.sd, ymax=frac_pos+frac_pos.sd),
         size=0.3, color=color_raw
      ) +
      geom_point(aes(fill=ge_min_samples_per_bin), shape=21, color='black') +
      scale_fill_manual(
         values=c(color_raw, color_filt),
         name='Samples per bin', labels=paste0(c('<','>='),min.samples.per.bin)
      ) +

      ## Filtered curve
      geom_path(data=subset(pd, ge_min_samples_per_bin), size=0.4, color=color_filt) +

      ## Axes
      scale_y_continuous(breaks=seq(0,1,0.2), name='Frac. target class samples in bin') +
      scale_x_continuous(breaks=seq(0,1,0.2), name='Prob. bin midpoint') +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +

      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         panel.grid.minor.y=element_blank(),
         legend.position='bottom'
      )

   if(show.labels){
      p <- p + geom_text(aes(label=label, y=label.ypos, hjust=label.hjust), angle=90, size=2.7)
   }

   return(p)
}
