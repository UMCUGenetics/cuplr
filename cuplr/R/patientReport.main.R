#' Plot patient report
#'
#' @description The patient report shows the predicted probabilities of each cancer type on a radar
#' plot, as well as the features contributing to the top predicted cancer types as bar plots
#'
#' @param probs A matrix where rows are samples, columns are cancer type prediction classes, and
#' cells are probabilities
#' @param feat.contrib A long form dataframe containing feature contribution data. This dataframe
#' contains the columns: sample, binary_rf, feature, contrib, feature_supplement
#' @param pred.report The input to this function can alternatively be a list containing the objects:
#' probs_adjusted, feat_contrib. This corresponds to `probs` and `feat.contrib`
#' @param sample.name Name of the sample in `pred.report` for which to make the patient report.
#' Required when input data is provided with `pred.report`
#' @param plot.title Plot title
#' @param top.n.pred.classes The number of panels to display that show the details for the top
#' predicted classes (including the exact probabilities, and feature contributions).
#' @param top.n.features The number features to show in the feature contribution plots in the
#' details pnaels
#' @param low.prob.thres If the top cancer type probability is lower than this value, the predicted
#' cancer type is considered to be uncertain. Multiple details panels are then shown.
#' Otherwise, only the panel for top predicted cancer type is shown
#' @param max.prob.diff Related to `low.prob.thres`. The lagged differences in probabilities for
#' the top cancer types (i.e. 1st to 2nd, 2nd to 3rd, 3rd to 4th, etc) is calculated. If e.g. the
#' differences between the 2nd to 3rd top cancer types is >=`max.prob.diff`, then details panels are
#' only shown for the 1st and 2nd top cancer types
#' @param rel.heights A numeric vector of length 2 specifying the relative heights of the radar plot
#' and feature contribution plots
#'
#' @return A ggplot grob
#' @export
#'
patientReport <- function(
   probs=NULL, feat.contrib=NULL,
   pred.report=NULL,
   sample.name=NULL, plot.title=sample.name,
   top.n.pred.classes=3, top.n.features=5,
   low.prob.thres=0.3, max.prob.diff=0.15,
   rel.heights=c(1,1)
){

   if(F){
      sample.name <- 'XXXXXXXX' ## unclear, gastric
      sample.name <- 'XXXXXXXX' ## Breast

      probs=pred_report.val$probs_adjusted[sample.name,]
      feat.contrib=subset(pred_report.val$feat_contrib, sample==sample.name)
      pred.report=NULL

      top.n.pred.classes=3
      top.n.features=5
      low.prob.thres=0.3
      max.prob.diff=0.15

      plot.title=''
   }

   ## Init --------------------------------
   require(ggplot2)

   if(!is.null(pred.report)){
      if(is.null(sample.name)){ stop('When input data is provided with `pred.report`, `sample.name` must also be provided') }
      probs <- pred.report$probs_adjusted[sample.name,]
      feat.contrib <- subset(pred.report$feat_contrib, sample==sample.name)
   } else if(is.null(probs) & is.null(feat.contrib)){
      stop('Input data must be provided to i) `probs` and `feat.contrib`, or ii) `pred.report`')
   }

   ## Colors --------------------------------
   ## Pred rank highlighting
   colors_dark <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999") ## Set1
   colors_light <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2") ## Pastel1

   #colors_dark <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
   #colors_light <- c("#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC")

   ##
   # feat_types <- unique(sub('[.].*$','',levels(feat.contrib$feature)))
   # feat.contrib$feature_type <- factor(
   #    sub('[.].*$','',feat.contrib$feature),
   #    feat_types
   # )
   #
   # colors_feat <- c(
   #    "#8DD3C7","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F", ##Excl: "#D9D9D9","#FFFFB3"
   #    "#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2",
   #    "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC"
   # )
   #
   # colors_feat <- structure(
   #    colors_feat[1:length(levels(feat.contrib$feature_type))],
   #    names=levels(feat.contrib$feature_type)
   # )

   ## Select classes --------------------------------
   probs_top <- probs[order(-probs)][1:top.n.pred.classes]

   if(probs_top[1]>=low.prob.thres){
      is_diff <- abs(c(0,diff(probs_top))) <= max.prob.diff ## Calc diff between the sorted probabilities
      probs_top <- probs_top[1:rle(is_diff)$lengths[1]] ## Select the top n classes where diff is still less than e.g. 0.15
   }

   ## Probs --------------------------------
   p_probs <- (function(){
      pd_probs <- as.data.frame(t(probs))
      pd_probs <- cbind(group=1, pd_probs)

      sel_colors <- structure(
         colors_dark[ 1:length(names(probs_top)) ],
         names=names(probs_top)
      )

      class_colors <- structure(
         rep('black',length(probs)),
         names=names(probs)
      )
      class_colors[names(sel_colors)] <- sel_colors

      ggradar(
         pd_probs,
         group.point.size=1.5, group.line.width=0.75,
         axis.label.size=3.5, axis.label.offset=1.1,
         axis.label.colors=class_colors,
         grid.label.size=5, gridline.label.offset=0,
         values.radar=c('0.0','0.5','1.0')
      ) +
         ggtitle(plot.title) +
         theme(
            plot.title=element_text(size=14)
         )
   })()

   ### Bar plot
   # pd_probs <- data.frame(class=names(probs), prob=unname(probs))
   # pd_probs <- pd_probs[order(-pd_probs$prob),]
   #
   # pd_probs$class_2 <- factor(pd_probs$class, rev(unique(pd_probs$class)))
   #
   # class_colors <- structure(
   #    rep('lightgrey',nrow(pd_probs)),
   #    names=pd_probs$class
   # )
   # class_colors[1:top.n.pred.classes] <- colors_dark[1:top.n.pred.classes]
   #
   # p_probs <- ggplot(pd_probs, aes(x=class_2, y=prob)) +
   #    geom_bar(aes(fill=class_2), stat='identity', color='black', size=0.3, show.legend=F) +
   #    geom_text(aes(y=prob+0.01, label=prob), hjust=0, size=3) +
   #    scale_fill_manual(values=class_colors) +
   #    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0,1,0.25), expand=c(0,0), name='Cancer type probability') +
   #
   #    coord_flip() +
   #    theme_bw() +
   #    theme(
   #       panel.grid.major.y=element_blank(),
   #       panel.grid.minor.y=element_blank(),
   #       panel.grid.minor.x=element_blank(),
   #       axis.title.y=element_blank()
   #    )

   ## Polar bar plot
   # pd_prob <- data.frame(
   #    pred_class=names(probs),
   #    prob=unname(probs)
   # )
   #
   # n_classes <- nrow(pd_prob)
   # pd_prob$index <- 1:n_classes
   #
   # pd_prob$pred_class <- factor(pd_prob$pred_class, unique(pd_prob$pred_class))
   #
   # ## Labels
   # pd_prob$angle <- 0:(n_classes-1) * 360/n_classes
   # pd_prob$label <- pd_prob$pred_class
   # pd_prob$label_hjust <- 0
   # pd_prob <- within(pd_prob,{
   #    label_hjust[round(angle) %in% c(0,180)] <- 0.5
   #    label_hjust[angle>180] <- 1
   # })
   #
   # plot_rotate_rad <- 0.5 * (360/n_classes) * (pi/180)
   #
   # coord_radar <- function (theta = "x", start = 0, direction = 1) {
   #    theta <- match.arg(theta, c("x", "y"))
   #    r <- if (theta == "x") "y" else "x"
   #    ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start,
   #            direction = sign(direction),
   #            is_linear = function(coord) TRUE)
   # }
   #
   # p_probs <- ggplot(pd_prob, aes(x=pred_class, y=prob)) +
   #
   #    ## Draw panel grid
   #    geom_segment(aes(x=pred_class, xend=pred_class, y=0, yend=1), color='lightgrey') +
   #    geom_hline(data=data.frame(), aes(yintercept=seq(0,1,0.25)), color='lightgrey') +
   #
   #    geom_segment(aes(x=pred_class, xend=pred_class, y=0, yend=prob), color='indianred', size=3) +
   #    geom_text(data=data.frame(x=1, y=seq(0,1,0.25)), aes(x=x, y=y, label=y), vjust=1.1) +
   #
   #    #geom_text(aes(y=label_pos, label=label, hjust=label_hjust)) +
   #    #geom_text(aes(y=1.1, label=pred_class, hjust=label_hjust)) +
   #    scale_y_continuous(expand=c(0, 0.2, 0, 0), limits=c(0, 1)) +
   #
   #    #coord_flip() +
   #    coord_polar(start = -plot_rotate_rad) +
   #    #coord_radar(start = -plot_rotate_rad) +
   #
   #    theme_bw() +
   #
   #    theme(
   #       panel.grid=element_blank(),
   #       panel.border=element_blank(),
   #       axis.text.x=element_text(size=12),
   #       axis.title.x=element_blank(),
   #       axis.title.y=element_blank(),
   #       axis.text.y=element_blank(),
   #       axis.ticks.y=element_blank()
   #    )


   ## Feat contrib --------------------------------
   feat_plots_xlims <- c(1-0.6, top.n.features+0.6)

   pd_feat <- feat.contrib[feat.contrib$binary_rf %in% names(probs_top),]

   ## Order by top predicted class and decreasing feat contrib
   pd_feat$binary_rf <- factor(pd_feat$binary_rf, names(probs_top))
   pd_feat <- pd_feat[order(pd_feat$binary_rf, -pd_feat$contrib),]

   ## Select top features
   pd_feat$feature_rank <- unlist(lapply(split(1:nrow(pd_feat), pd_feat$binary_rf),function(i){
      1:length(i)
   }), use.names=F)
   pd_feat <- pd_feat[pd_feat$feature_rank<=top.n.features,]

   p_contrib <- ggplot(pd_feat, aes(x=-feature_rank, y=contrib)) +
      coord_flip() +
      facet_grid(binary_rf~., scales='fixed', switch='y') +

      geom_bar(stat='identity', fill='lightgrey', color='darkgrey', size=0.3, width=0.8) +
      geom_text(aes(y=0, label=feature), size=3.5, angle=0, hjust=0) +
      ylab('Feature contribution') +
      scale_x_continuous(expand=c(0, 0), limits=-rev(feat_plots_xlims)) +

      theme_bw() +
      theme(
         panel.grid.major.y=element_blank(),
         #panel.grid.minor.y=element_blank(),
         panel.grid.minor.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()
      )


   ## Edit strip colors
   colorizeStrips <- function(p){
      #p=p_contrib
      g <- ggplot_gtable(ggplot_build(p))
      strip_both <- which(grepl('strip-', g$layout$name))
      fills <- structure(
         colors_light[ 1:length(names(probs_top)) ],
         names=names(probs_top)
      )

      counter <- 1
      for (i in strip_both) {
         j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
         g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[counter]
         counter <- counter+1
      }
      return(g)
   }

   p_contrib <- colorizeStrips(p_contrib)

   ## Cohort summary stats --------------------------------
   feat_meta_colnames <- c('binary_rf','feature','feature_rank')
   feat_stat_colnames <- c('avg_case','avg_ctrl','value')

   pd_feat_avg <- reshape2::melt(
      pd_feat[,c(feat_meta_colnames, feat_stat_colnames)],
      id.vars=feat_meta_colnames
   )

   ## Center data based on the avg values of the ctrl cohort
   norm_stats <- t(apply(
      pd_feat[,feat_stat_colnames],
      1,
      function(x){ (x-min(x))/(max(x)-min(x)) }
   ))
   norm_stats <- as.data.frame(norm_stats)
   shift_factor <- 1 - norm_stats$avg_ctrl
   norm_stats <- norm_stats + shift_factor

   colnames(norm_stats) <- paste0('norm.',colnames(norm_stats))
   norm_stats <- cbind(pd_feat[,feat_meta_colnames], norm_stats)
   norm_stats <- reshape2::melt(norm_stats,id.vars=feat_meta_colnames)

   pd_feat_avg$value_norm <- norm_stats$value
   rm(norm_stats)

   ## Stagger lollipops
   dodge_width <- 0.2
   pd_feat_avg$xpos <- pd_feat_avg$feature_rank
   pd_feat_avg$xpos[pd_feat_avg$variable=='avg_case'] <- pd_feat_avg$xpos[pd_feat_avg$variable=='avg_case'] + dodge_width
   pd_feat_avg$xpos[pd_feat_avg$variable=='value'] <- pd_feat_avg$xpos[pd_feat_avg$variable=='value'] - dodge_width

   ## Reverse x values for coord_flip
   pd_feat_avg$xpos <- -pd_feat_avg$xpos

   ## Calc label position
   label_offset <- 0.05
   pd_feat_avg$label_ypos <- pd_feat_avg$value_norm
   pd_feat_avg$label_hjust <- 1

   pd_feat_avg <- within(pd_feat_avg,{
      label_ypos[label_ypos >1] <- label_ypos[label_ypos >1] + label_offset
      label_ypos[label_ypos<=1] <- label_ypos[label_ypos<=1] - label_offset
      label_hjust[label_ypos>1] <- 0
   })

   ## Reorder value type
   pd_feat_avg$variable <- factor(pd_feat_avg$variable, c('value','avg_case','avg_ctrl'))
   value_type_names <- c('Sample','Mean, cancer type','Mean, other')

   ## Format label
   pd_feat_avg$value <- signif(pd_feat_avg$value,3)
   #log10(abs(pd_feat_avg$value))
   #as.character(pd_feat_avg$value)

   p_feat_avg <- ggplot(pd_feat_avg, aes(x=xpos, y=value_norm)) +
      facet_grid(binary_rf~., scales='fixed') +

      ## Baseline
      geom_segment(
         data=subset(pd_feat_avg, variable=='avg_ctrl'),
         mapping=aes(x=xpos-dodge_width-0.1, xend=xpos+dodge_width+0.1, y=1, yend=1),
         show.legend=F
      ) +

      ## Data points
      geom_segment(aes(x=xpos, xend=xpos, y=1, yend=value_norm, color=variable), show.legend=F) +
      geom_point(aes(shape=variable, size=variable, color=variable), show.legend=T) +
      geom_text(aes(y=label_ypos, label=value, hjust=label_hjust, color=variable), show.legend=F, size=2.5) +

      ## Legend
      scale_color_manual(values=c(avg_ctrl='black',avg_case='darkgrey',value='#F58225'), labels=value_type_names) +
      scale_shape_manual(values=c(avg_ctrl=124,avg_case=16,value=16), labels=value_type_names) +
      scale_size_manual(values=c(avg_ctrl=0,avg_case=2,value=2), labels=value_type_names) +
      guides(shape=F, size=F) +

      ## Axes
      scale_y_continuous(expand=c(0.1, 0.1), name='Feature value (variable units)') +
      scale_x_continuous(expand=c(0, 0), limits=-rev(feat_plots_xlims)) +

      coord_flip() +
      theme_bw() +
      theme(
         panel.grid.major.y=element_blank(),
         #panel.grid.minor.y=element_line(color='grey'),
         panel.grid.major.x=element_blank(),
         panel.grid.minor.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.position='bottom',
         legend.title=element_blank(),
         legend.direction='horizontal',
         strip.background = element_blank(),
         strip.text = element_blank()
      )

   #p_feat_avg <- colorizeStrips(p_feat_avg)

   ## Merge plots --------------------------------
   cowplot::plot_grid(
      p_probs,
      cowplot::plot_grid(
         p_contrib, p_feat_avg,
         nrow=1, align='h', axis='tblr'
      ),
      ncol=1,
      rel_heights=rel.heights
   )


}
