#' Calculate cumulative accuracy and optimum probability cutoffs
#'
#' @description Calculate cumulative accuracy, and from this calculute optimum probability cutoffs
#' for each binary random forest. The output of this function also serves as the plot data from
#' `plotSortedProbs()`
#'
#' @param actual A vector of the actual classes
#' @param predicted A vector of the predicted classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param min.frac.correct For each (actual) class, the threshold for the minimum fraction of
#' samples correctly predicted. The probability at which this fraction occurs is assigned as the
#' optimum probability cutoff
#'
#' @return A dataframe where samples are ordered by their actual class (alphabetical), followed by
#' the predicted probability for the actual class (descending). The `is_opt_thres` column indicates
#' the rows corresponding to the probabilities with the optimum thresholds
#' @export
#'
calcCumAcc <- function(actual, predicted, probs, min.frac.correct=1){
   if(F){
      devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CUPs_classifier/processed/cuplr/commonUtils/')
      actual=pred_reports.cv$HMF$responses_actual
      predicted=pred_reports.cv$HMF$responses_pred
      probs=pred_reports.cv$HMF$probs_adjusted
   }

   ## Init --------------------------------
   if(length(unique(c(length(actual),length(predicted),nrow(probs))))!=1){
      stop('length(actual), length(predicted), and nrow(probs) must be equal')
   }

   df <- data.frame(actual,predicted,row.names=NULL)
   df <- cbind(sample=1:nrow(df), df)

   ## Get prob of target cancer type --------------------------------
   probs_melt <- probs; rownames(probs_melt) <- NULL
   probs_melt <- reshape2::melt(probs_melt)
   colnames(probs_melt) <- c('sample','binary_rf','prob')
   df$prob_target_class <- probs_melt$prob[
      match(
         paste0(df$sample, as.character(df$actual)),
         paste0(probs_melt$sample, probs_melt$binary_rf)
      )
   ]
   rm(probs_melt)

   ## Calculate cum. no. of correct samples --------------------------------
   df <- df[order(df$actual, -df$prob_target_class),]
   df$sample <- factor(df$sample, unique(df$sample))

   ## Total no. of samples per cancer type
   class_counts <- table(df$actual)
   df$n_in_class <- as.integer(class_counts[as.character(df$actual)])
   rm(class_counts)

   ## Cum. no. of samples correct
   df$is_correct_pred <- as.character(df$actual)==as.character(df$predicted)
   l <- split(df$is_correct_pred, df$actual)
   df$n_correct_pred_cum <- unlist(lapply(l, function(i){
      #i=l[[1]]
      n_correct <- rep(length(i), length(i))
      n_wrong_cum <- cumsum(!i)
      n_correct - n_wrong_cum
   }))
   df$is_correct_pred <- NULL; rm(l)

   df$frac_correct_pred_cum <- df$n_correct_pred_cum / df$n_in_class

   ## Determine optimal cutoff --------------------------------
   l <- split(df$frac_correct_pred_cum, df$actual)
   df$is_opt_thres <- unlist(lapply(l, function(i){
      #i=l[[4]]

      is_opt_thres <- rep(FALSE, length(i))
      le_min_frac_correct <- i>=min.frac.correct
      if(all(!le_min_frac_correct)){
         is_opt_thres[ max(which(i==max(i))) ] <- TRUE
         warning('Some classes did not have `min.frac.correct`>=',min.frac.correct)
      } else {
         is_opt_thres[ max(which(le_min_frac_correct)) ] <- TRUE
      }

      return(is_opt_thres)
   }))
   rm(l)

   #subset(df, actual=='HeadAndNeck_ACC')

   ## Get sample rank for plotting --------------------------------
   #df <- df[order(df$actual, df$prob_target_class),]

   l <- split(df$sample, df$actual)
   df$rank <- unlist(lapply(l, function(i){
      0:(length(i)-1)
   }))
   rm(l)

   df$sample <- factor(df$sample, unique(df$sample))

   rownames(df) <- NULL
   return(df)
}

####################################################################################################
#' Plot sort probabilities
#'
#' @param actual A vector of the actual classes
#' @param predicted A vector of the predicted classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param df Instead of providing `actual`,`predicted`, and `probs`, the output from
#' `calcClassifThres()` can be provided to `df`
#' @param min.frac.correct For each (actual) class, the threshold for the minimum fraction of
#' @param show.opt.thres Show the lines/labels at the optimum probability cutoff for each actual
#' class?
#' @param label.vjust A numeric value indicating the vertical justification of the `show.opt.thres`
#' label. Alternative, if 'auto', the vjust will automatically be calculated.
#' @param facet.ncol Number of facet columns. If 'auto', this value will be
#' `round(sqrt(n_classes))`
#'
#' @description For each actual class, make barplots of sort sample prediction probabilities
#'
#' @return A ggplot object
#' @export
#'
plotSortedProbs <- function(
   actual=NULL, predicted=NULL, probs=NULL, df=NULL, min.frac.correct=1,
   show.opt.thres=T, label.vjust='auto', facet.ncol='auto'
){

   if(F){
      actual=pred_report.cv$responses_actual
      predicted=pred_report.cv$responses_pred
      probs=pred_report.cv$probs_adjusted
      min.frac.correct=1
      show.opt.thres=T
      label.vjust='auto'
   }

   ## Init --------------------------------
   require(ggplot2)

   if(is.null(df)){
      df <- calcCumAcc(
         actual=actual, predicted=predicted, probs=probs,
         min.frac.correct=min.frac.correct
      )
   }
   #df <- df[df$actual %in% c('Bile_Gallbladder','Breast','Lung_NSC','Glioma','Lymphoid'),]

   ## Predicted class colors --------------------------------
   color_pal <- c(
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", ## Set3
      "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", ## Set2
      "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", ## Pastel1
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999" ## Set1
   )

   uniq_classes <- sort(unique(as.character(df$predicted)))
   class_colors <- structure(
      color_pal[1:length(uniq_classes)],
      names=uniq_classes
   )

   ## Misc plot data --------------------------------
   ## Optimum cutoff line
   thresholds <- df[df$is_opt_thres,]

   if(label.vjust=='auto'){
      thresholds$label_vjust <- 0
      thresholds$label_vjust[thresholds$prob_target_class>=0.25] <- 1
   } else {
      thresholds$label_vjust <- label.vjust
   }

   ## Close off geom step
   l <- split(df$rank, df$actual)
   l <- l[sapply(l, length)!=0]
   df$is_last_rank <- unlist(lapply(l, function(i){
      i==max(i)
   }), use.names=F)
   last_rank <- df[df$is_last_rank,]

   ## Plot --------------------------------
   if(facet.ncol=='auto'){
      facet.ncol <- round(sqrt(length(unique(actual))))
   }

   p <- ggplot(df, aes(x=rank, y=prob_target_class)) +
      facet_wrap(~actual, scales='free_x', ncol=facet.ncol) +

      ## Bars + outline
      geom_bar(aes(fill=predicted), stat='identity', width=1) +

      geom_step(aes(x=rank-0.5), size=0.3) +
      geom_segment(
         data=last_rank,
         aes(x=rank-0.5, xend=rank+0.5, y=prob_target_class, yend=prob_target_class),
         size=0.3
      ) +
      scale_fill_manual(name='Predicted class', values=class_colors) +
      scale_y_continuous(name='Prob. of target class', limits=c(0,1)) +
      scale_x_continuous(name='Rank') +

      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         panel.grid.minor.y=element_blank(),
         axis.ticks.x=element_blank()
      )

   ## Indicate optimum threshold lines
   if(show.opt.thres){
      p <- p +
         geom_hline(
            data=thresholds, aes(yintercept=prob_target_class),
            size=0.3, linetype='dashed'
         ) +
         geom_label(
            data=thresholds, aes(x=0, y=prob_target_class, label=prob_target_class, vjust=label_vjust),
            hjust=0, size=3.5, alpha=0.5
         )
   }

   return(p)
}
