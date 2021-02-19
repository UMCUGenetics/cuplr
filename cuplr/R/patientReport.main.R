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
#' @param feat.contrib.scales Should scales be fixed ("fixed", the default), free ("free"), or free
#' in one dimension ("free_x", "free_y")?
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
   feat.contrib.scales='fixed'
){

   if(F){
      #subset(metadata$training, training_group!=1)
      sample.name <- 'XXXXXXXX' ## Prostate

      probs=pred_report$probs_adjusted[sample.name,]
      feat.contrib=subset(pred_report$feat_contrib, sample==sample.name)

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

   ## Radar plot --------------------------------
   p_radar <- (function(){
      pd <- as.data.frame(t(probs))
      pd <- cbind(group=1, pd)

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
         pd,
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

   ## Feat contrib --------------------------------
   p_feat_contrib <- (function(){

      pd <- feat.contrib[feat.contrib$binary_rf %in% names(probs_top),]

      ## Order by top predicted class and decreasing feat contrib
      pd$binary_rf <- factor(pd$binary_rf, names(probs_top))
      pd <- pd[order(pd$binary_rf, -pd$contrib),]

      ## Select top features
      pd$index <- unlist(lapply(split(1:nrow(pd), pd$binary_rf),function(i){
         1:length(i)
      }), use.names=F)
      pd <- pd[pd$index<=top.n.features,]

      ## Add top binary RF probs
      pd$prob <- probs[ match(pd$binary_rf, names(probs)) ]
      pd$strip_title <- paste0(pd$binary_rf,' (Probability=',pd$prob,')')
      pd$strip_title <- factor(pd$strip_title, unique(pd$strip_title))

      ## Add feature supplement to name
      pd$feature_supplement <- paste0(' (',pd$feature_supplement,')')
      pd$feature_supplement[pd$feature_supplement==' ()'] <- ''
      pd$label <- paste0(as.character(pd$feature),pd$feature_supplement)

      ## Main
      p <- ggplot(pd, aes(x=-index, y=contrib)) +
         coord_flip() +
         facet_wrap(~strip_title, ncol=1, scales=feat.contrib.scales) +

         geom_bar(stat='identity', fill='lightgrey', color='darkgrey', size=0.3) +
         geom_text(aes(y=0, label=label), size=3.5, angle=0, hjust=0) +
         ylab('Feature contribution') +

         theme_bw() +
         theme(
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
         )

      ## Edit strip colors
      g <- ggplot_gtable(ggplot_build(p))
      strip_both <- which(grepl('strip-', g$layout$name))
      fills <- rev(structure(
         colors_light[ 1:length(names(probs_top)) ],
         names=names(probs_top)
      ))

      counter <- 1
      for (i in strip_both) {
         j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
         g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[counter]
         counter <- counter+1
      }
      return(g)
   })()

   out <- gridExtra::arrangeGrob(
      p_radar,
      p_feat_contrib,
      nrow=1, widths=c(2, 1)
   )
   return(out)
}
