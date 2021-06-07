#' Plot patient report
#'
#' @description The patient report shows the predicted probabilities of each cancer type on a radar
#' plot, as well as the features contributing to the top predicted cancer types as bar plots
#'
#' @param probs A numeric vector where the names are the predicted classes
#' @param feat.contrib A long form dataframe containing feature contribution data. This dataframe
#' contains the columns: sample, binary_rf, feature, contrib, feature_supplement
#' @param sample.name Name of the sample in `report` for which to make the patient report.
#' Required when input data is provided with `report`
#' @param plot.title Plot title
#' @param top.n.pred.classes The number of panels to display that show the details for the top
#' predicted classes (including the exact probabilities, and feature contributions).
#' @param top.n.features The number features to show in the feature contribution plots in the
#' details panels
#' @param rel.widths A numeric vector of length 2 specifying the relative widths of the probability
#' plot and feature contribution plots
#' @param prob.thres.min The feature contribution panels will be shown for predicted classes with a
#' probability higher than this threshold
#' @param prob.thres.rel.diff The feature contribution panels will be shown for the top predicted
#' classes, as long as the relative difference between the 1st and 2nd classes, 2nd and 3rd classes,
#' etc, are higher than this threshold. Note that `prob.thres.min` overrides this threshold
#'
#' @return A ggplot grob
#' @export
#'
patientReport <- function(
   probs=NULL, feat.contrib=NULL,
   sample.name=NULL, plot.title=sample.name,
   top.n.pred.classes=3, top.n.features=5,
   prob.thres.min=0.1, prob.thres.rel.diff=0.4,
   rel.widths=c(2,1)
){

   # if(F){
   #    probs=pred_report$prob_scaled
   #    feat.contrib=pred_report$feat_contrib
   #    sample.name <- 'DO36107' ## Uncertain
   #    plot.title=sample.name
   #
   #    top.n.pred.classes=3
   #    top.n.features=5
   #    prob.thres.min=0.1
   #    prob.thres.rel.diff=0.4
   #    rel.widths=c(2,1)
   # }

   ## Init --------------------------------
   require(ggplot2)

   if(is.null(rownames(probs))){ rownames(probs) <- 1:nrow(probs) }
   if(is.null(sample.name) & nrow(probs)==1){
      sample.name <- rownames(probs)
      probs <- probs[1,]
   }
   if(is.null(sample.name) & nrow(probs)>=1){
      stop('`probs` contains multiple samples. Please specify a `sample.name`')
   }
   probs <- probs[sample.name,]

   if(!is.numeric(probs) | !is.vector(probs) | length(probs)==0){
      stop('`probs` must be a named numeric vector')
   }

   if(!is.data.frame(feat.contrib)){
      stop('`feat.contrib` must be a long form dataframe')
   }

   ## Colors --------------------------------
   ## Pred rank highlighting
   colors_dark <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999") ## Set1
   colors_light <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2") ## Pastel1

   ## Select classes --------------------------------
   probs_top <- probs[order(-probs)][1:top.n.pred.classes]

   lagged_rel_diff <-
      probs_top /
      c(
         probs_top[1],
         probs_top[ (1:length(probs_top)-1) ]
      )
   lagged_rel_diff[is.na(lagged_rel_diff)] <- 0

   ## Use trycatch to prevent -Inf when also probs don't meet the thresholds
   last_sel_class.rel_diff <- tryCatch(
      expr = { max(which(lagged_rel_diff>=prob.thres.rel.diff)) },
      warning = function(w){ return(1) }
   )

   last_sel_class.min_prob <- tryCatch(
      expr = { max(which(probs_top>=prob.thres.min)) },
      warning = function(w){ return(1) }
   )

   last_sel_class <- max(last_sel_class.rel_diff, last_sel_class.min_prob)

   probs_top <- probs_top[ 1:last_sel_class ]

   ## Probs --------------------------------
   pd_probs <- data.frame(class=names(probs), prob=unname(probs))
   pd_probs <- pd_probs[order(-pd_probs$prob),]
   pd_probs$class <- factor(pd_probs$class, pd_probs$class)

   pd_probs$label <- format(round(pd_probs$prob, 3), nsmall=3)
   pd_probs$label[pd_probs$prob==0] <- '0'

   ## Colors
   ## Recycle colors when many classes need to be colored
   color_indexes <- rep(
      1:length(colors_light),
      ceiling(nrow(pd_probs)/length(colors_light))
   )
   color_indexes <- color_indexes[1:nrow(pd_probs)]

   line_color <- structure(
      colors_dark[color_indexes],
      names=as.character(pd_probs$class)
   )
   line_color[pd_probs$prob==0] <- 'lightgrey'

   label_fill <- structure(
      colors_light[color_indexes],
      names=as.character(pd_probs$class)
   )
   label_fill[pd_probs$prob==0] <- 'lightgrey'

   ## Main
   p_probs <- ggplot(pd_probs, aes(x=class, y=prob)) +

      geom_segment(aes(x=class, xend=class, y=0, yend=prob, color=class), size=2, show.legend=F) +
      scale_color_manual(values=line_color) +

      geom_label(aes(fill=class, label=label), size=3, show.legend=F, label.padding=unit(0.15,'lines')) +
      scale_fill_manual(values=label_fill) +

      scale_y_continuous(name='Cancer type probability', limits=c(0, 1), breaks=seq(0,1,0.2), expand=c(0.05, 0.05)) +
      scale_x_discrete(name='Cancer type', limits=rev) +

      coord_flip() +
      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_text(size=10),
         axis.text.x=element_text(size=10)
      )

   ## Feat contrib --------------------------------
   pd_feat <- feat.contrib[feat.contrib$binary_rf %in% names(probs_top),]

   if(!all(names(probs_top) %in% unique(pd_feat$binary_rf))){
      stop('`feat.contrib` is missing data from the top predicted classes')
   }

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
      facet_wrap(binary_rf~., scales='fixed', ncol=1) +

      geom_bar(stat='identity', fill='lightgrey', color='darkgrey', size=0.3, width=0.8) +
      geom_text(aes(label=feature), y=max(pd_feat$contrib)*0.01, size=3.5, angle=0, hjust=0) +
      ylab('Feature contribution') +

      theme_bw() +
      theme(
         panel.grid.major.y=element_blank(),
         panel.grid.minor.x=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank()
      )

   ## Edit strip colors
   colorizeStrips <- function(p, colors){

      # if(F){
      #    p=p_contrib
      #
      #    colors=structure(
      #       colors_light[ 1:length(names(probs_top)) ],
      #       names=names(probs_top)
      #    )
      # }

      g <- ggplot_gtable(ggplot_build(p))
      names(g$grobs) <- g$layout$name
      strip_ids <- grep('strip-', g$layout$name, value=T)
      strip_ids <- strip_ids[
         order(sapply(strsplit(strip_ids,'-'), function(i){ i[length(i)] }))
      ]
      strip_ids <- na.exclude(strip_ids[1:length(colors)])

      counter <- 0
      for(i in strip_ids){
         #i=strip_ids[[1]]
         counter <- counter + 1
         object <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
         g$grobs[[i]]$grobs[[1]]$children[[object]]$gp$fill <- colors[counter]
      }

      p <- ggplotify::as.ggplot(g)
      return(p)
   }

   p_contrib <- colorizeStrips(
      p_contrib,
      structure(
         colors_light[ 1:length(names(probs_top)) ],
         names=names(probs_top)
      )
   )

   ## Combine plots --------------------------------
   plot_title <- NULL
   if(!is.null(plot.title)){
      plot_title <- grid::textGrob(plot.title, x=0.02, hjust=0)
   }

   gridExtra::arrangeGrob(
      p_probs, p_contrib, nrow=1,
      top=plot_title, widths=rel.widths
   )
}
