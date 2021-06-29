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
#' @param prob.thres.min The feature contribution panels will be shown for predicted classes with a
#' probability higher than this threshold
#' @param prob.thres.rel.diff The feature contribution panels will be shown for the top predicted
#' classes, as long as the relative difference between the 1st and 2nd classes, 2nd and 3rd classes,
#' etc, are higher than this threshold. Note that `prob.thres.min` overrides this threshold
#' @param which.plots A character vector indicating which plots to show. Can be 'probs',
#' 'feat.values', and/or 'feat.contribs'
#' @param rel.widths A numeric vector of the same length as `which.plots` specifying the relative
#' widths of the subplots
#'
#' @return A cowplot of ggplot objects
#' @export
#'
patientReport <- function(
   probs=NULL, feat.contrib=NULL,
   sample.name=NULL, plot.title=sample.name,
   top.n.pred.classes=3, top.n.features=5,
   prob.thres.min=0.1, prob.thres.rel.diff=0.4,
   which.plots=c('probs','feat.values'), rel.widths=c(1.4, 1)
){

   if(F){
      ##
      probs=pred_reports$holdout$prob_scaled
      feat.contrib=pred_reports$holdout$feat_contrib
      sample.name='DO36021' ## 2 top classes

      ##
      probs=pred_reports$CV$prob_scaled
      feat.contrib=pred_reports$CV$feat_contrib
      sample.name='DO52615' ## 3 top classes

      plot.title=sample.name

      top.n.pred.classes=3
      top.n.features=5
      prob.thres.min=0.1
      prob.thres.rel.diff=0.4
      rel.widths=c(2,1)
   }

   ## Init ================================
   require(ggplot2)

   ##
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

   ##
   if(!is.data.frame(feat.contrib)){
      stop('`feat.contrib` must be a long form dataframe')
   }
   feat.contrib <- subset(feat.contrib, sample==sample.name)

   ## Colors ================================
   ## Pred rank highlighting
   colors_dark  <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999") ## Set1
   colors_light <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2") ## Pastel1

   ## Select classes ================================
   probs_top <- probs[order(-probs)][1:top.n.pred.classes]

   lagged_rel_diff <-
      probs_top /
      c(
         probs_top[1],
         probs_top[ (1:length(probs_top)-1) ]
      )
   lagged_rel_diff[is.na(lagged_rel_diff)] <- 0

   ## Use trycatch to prevent -Inf when probs don't meet the thresholds
   ## When -Inf occurs, always return the top class prob
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


   ## Probs ================================
   ## Put probs into long form dataframe --------------------------------
   pd_probs <- data.frame(class=names(probs), prob=unname(probs))
   pd_probs <- pd_probs[order(-pd_probs$prob),]
   pd_probs$class <- factor(pd_probs$class, pd_probs$class)

   pd_probs$prob_round <- round(pd_probs$prob, 3)
   pd_probs$label <- format(pd_probs$prob_round, nsmall=3)
   pd_probs$label[pd_probs$prob_round==0] <- '0'

   ## Colors --------------------------------
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
   line_color[pd_probs$prob_round==0] <- 'lightgrey'

   label_fill <- structure(
      colors_light[color_indexes],
      names=as.character(pd_probs$class)
   )
   label_fill[pd_probs$prob_round==0] <- 'lightgrey'

   ## Main --------------------------------
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

   if(!is.null(plot.title)){
      p_probs <- p_probs + ggtitle(plot.title)
   }


   ## Feat contrib ================================
   ## Prep data --------------------------------
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

   ## Assign a unique feature label for each class
   ## Ensures feature rank is preserved for each class
   pd_feat$feature_label <- paste0(pd_feat$binary_rf, '::', pd_feat$feature)
   pd_feat$feature_label <- factor(pd_feat$feature_label, rev(unique(pd_feat$feature_label)))

   p_contrib <- ggplot(pd_feat, aes(x=feature_label, y=contrib)) +
      coord_flip() +
      facet_wrap(binary_rf~., scales='free_y', ncol=1) +

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

   ## Add blank title to make sure plots are aligned
   if(!is.null(plot.title)){
      p_contrib <- p_contrib + ggtitle(' ')
   }

   ## Edit strip colors --------------------------------
   colorizeStrips <- function(p, colors){
      # if(F){
      #    p=p_contrib
      #    colors=structure( colors_light[ 1:length(names(probs_top)) ], names=names(probs_top) )
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

      return(g)

      ## Using ggplotify messes up alignment with cowplot::plot_grid()
      #p <- ggplotify::as.ggplot(g)
      #return(p)
   }

   strip_colors <- structure(
      colors_light[ 1:length(names(probs_top)) ],
      names=names(probs_top)
   )

   p_contrib <- colorizeStrips(p_contrib, strip_colors)


   ## Feat values ================================
   ## Prep data --------------------------------
   feat_value_types <- data.frame(
      matrix(c(
         'min_all',      'black',     108, 3,   1,
         'max_all',      'black',     108, 3,   1,
         'avg_ctrl',     'black',     108, 3,   1,
         'avg_case',     'indianred', 16,  2, 1,
         'sample_value', 'indianred', 21,  2.5, 0
      ),
      dimnames=list(NULL, c('type','color','shape','size','alpha')),
      ncol=5, byrow=T)
   )

   constructAesValues <- function(colname){
      values <- feat_value_types[,colname]
      if(grepl('^\\d',values[1])){
         values <- as.numeric(values)
      }
      structure(
         values,
         names=feat_value_types$type
      )
   }

   ## Calculate min/max feature value for each feature
   pd_feat$stat_min <- apply(pd_feat[,feat_value_types$type],1,min)
   pd_feat$stat_max <- apply(pd_feat[,feat_value_types$type],1,max)

   ## Conver to long form dataframe
   pd_feat_values <- lapply(feat_value_types$type, function(i){
      cbind(
         pd_feat[,c('binary_rf','feature','feature_label','feature_rank','stat_min','stat_max','avg_metric')], ## Fixed columns
         value_type=i,
         value=pd_feat[,i] ## Variable (value) columns
      )
   })
   pd_feat_values <- do.call(rbind, pd_feat_values)

   ## Force order of plotting points
   pd_feat_values$value_type <- factor(pd_feat_values$value_type, feat_value_types$type)

   ## Log transform features with a wide range
   pd_feat_values$is_wide_range <- with(pd_feat_values,{
      stat_max-stat_min > 100 & avg_metric != 'prop'
   })
   pd_feat_values$value.trans <- pd_feat_values$value
   pd_feat_values$stat_min.trans <- pd_feat_values$stat_min
   pd_feat_values$stat_max.trans <- pd_feat_values$stat_max

   pd_feat_values <- within(pd_feat_values,{
      value.trans[is_wide_range] <- log10(value.trans[is_wide_range]+1)
      stat_min.trans[is_wide_range] <- log10(stat_min.trans[is_wide_range]+1)
      stat_max.trans[is_wide_range] <- log10(stat_max.trans[is_wide_range]+1)
   })

   ## Scale feature values to arbitrary units ranging from 0 to 1
   scale0to1 <- function(x, x.min, x.max){
      out <- (x-x.min)/(x.max-x.min)

      ## Clip out of bounds feature values to (0,1)
      out[out>1] <- 1
      out[out<0] <- 0

      return(out)
   }
   pd_feat_values$value.scaled <- with( pd_feat_values, scale0to1(value.trans, stat_min.trans, stat_max.trans) )

   ## Top to bottom feature order. 2 line spacing betwen features
   pd_feat_values$feature_rank <- -pd_feat_values$feature_rank * 2

   ## Add plotting values to data --------------------------------
   ## Make values human readable
   formatNumHuman <- function(x, avg.metric, sig.digits=3){
      x <- signif(x, sig.digits)
      out <- sapply(x, function(i){
         if(i<1e3){ return(i) }
         if(i<1e6){ return(paste0(i/1e3,'K')) }
         if(i<1e9){ return(paste0(i/1e6,'M')) }
         if(i<1e12){ return(paste0(i/1e9,'T')) }
         return(i)
      })

      out[avg.metric=='prop'] <- paste0( x[avg.metric=='prop']*100, '%' )

      return(out)
   }
   pd_feat_values$value_label <- with(pd_feat_values, formatNumHuman(value, avg_metric, sig.digits=2))

   ## Add hjust to ensure label doesnt go past min/max
   pd_feat_values$hjust <- 0.5
   pd_feat_values$hjust[pd_feat_values$value.scaled>=0.95] <- 1
   pd_feat_values$hjust[pd_feat_values$value.scaled<=0.05] <- 0

   ## Min/max value labels --------------------------------
   ## Get original min/max values
   pd_feat_values.range <- reshape2::dcast(
      data=subset(pd_feat_values, value_type %in% c('min_all','max_all')),
      formula=binary_rf+feature_rank+feature+is_wide_range+avg_metric~value_type,
      value.var='value'
   )
   pd_feat_values.range <- reshape2::melt(
      data=pd_feat_values.range,
      measure.vars=c('min_all','max_all')
   )

   pd_feat_values.range$value_label <- formatNumHuman(
      pd_feat_values.range$value,
      pd_feat_values.range$avg_metric
   )

   ## Assign justification and offsets
   pd_feat_values.range$value.scaled <- -0.02
   pd_feat_values.range$hjust <- 0
   pd_feat_values.range <- within(pd_feat_values.range,{
      value.scaled[variable=='max_all'] <- 1.02
      hjust[variable=='min_all'] <- 1
   })

   ## Segment between case and ctrl averages --------------------------------
   pd_feat_values.direction <- reshape2::dcast(
      data=subset(pd_feat_values, value_type %in% c('avg_case','avg_ctrl')),
      formula=binary_rf+feature_rank~value_type,
      value.var='value.scaled'
   )
   pd_feat_values.direction <- within(pd_feat_values.direction,{
      min <- pmin(avg_case, avg_ctrl)
      max <- pmax(avg_case, avg_ctrl)
   })

   ## Plot --------------------------------
   legend_labels <- c('Avg. in other cancer types','Avg. in target cancer type','Value in patient')

   p_feat_values <- ggplot(pd_feat_values, aes(x=feature_rank, y=value.scaled)) +
      coord_flip(ylim=c(-0.05, 1.1)) +
      facet_wrap(binary_rf~., scales='free_y', ncol=1) +

      ## Min/max segment
      geom_segment(
         data=pd_feat_values.range,
         mapping=aes(x=feature_rank, xend=feature_rank, y=0, yend=1),
         color='grey', size=0.2
      ) +

      ## Min/max labels
      geom_text(
         data=pd_feat_values.range,
         mapping=aes(x=feature_rank, y=value.scaled, label=value_label, hjust=hjust),
         vjust=0.4, size=3, color='grey60', fontface='bold'
      ) +

      ## Feature label
      geom_text(
         data=pd_feat_values.range,
         mapping=aes(x=feature_rank+0.8, y=0, label=feature),
         hjust=0, size=3, color='grey35'
      ) +

      ## Segment between case and ctrl averages
      geom_segment(
         data=pd_feat_values.direction,
         mapping=aes(x=feature_rank, xend=feature_rank, y=min, yend=max),
         color='indianred'
      ) +

      ## Feature stats points. Show value for patient.
      geom_point(
         data=subset(pd_feat_values, value_type %in% c('avg_ctrl','avg_case','sample_value')),
         aes(color=value_type, shape=value_type, size=value_type, alpha=value_type),
         fill='white'
      ) +
      geom_label(
         data=subset(pd_feat_values, value_type %in% c('sample_value')),
         aes(color=value_type, label=value_label, hjust=hjust),
         size=2.5, label.padding=unit(0.15,'lines'), show.legend=F
      ) +
      scale_color_manual(values=constructAesValues('color'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
      scale_shape_manual(values=constructAesValues('shape'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
      scale_size_manual (values=constructAesValues('size'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
      scale_alpha_manual(values=constructAesValues('alpha'), labels=legend_labels, guide=guide_legend(reverse=TRUE, override.aes=list(alpha=NA))) + ## alpha=NA: so that point on legend is drawn but point on plot not drawn

      ##
      scale_y_continuous(
         name='Feature value (variable units)',
         breaks=c(0, 0.25, 0.5, 0.75, 1),
         labels=c('Min','','','','Max')
      ) +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         legend.position='bottom',
         legend.direction='vertical',
         legend.margin=margin(0,0,0,0),
         legend.box.margin=margin(-10,-10, 0,-10),
         legend.key.height=unit(10,'pt'),
         legend.title=element_blank()
      )

   ## Add blank title to make sure plots are aligned
   if(!is.null(plot.title)){
      p_feat_values <- p_feat_values + ggtitle(' ') ## Add blank title to make sure plots are aligned
   }

   ## Add strip colors
   p_feat_values <- colorizeStrips(p_feat_values, strip_colors)


   ## Combine plots ================================
   #which.plots=c('probs','feat.contribs','feat.values')

   l_plots <- list()
   if('probs' %in% which.plots){
      l_plots$probs <- p_probs
   }
   if('feat.contribs' %in% which.plots){
      l_plots$feat.contribs <- p_contrib
   }
   if('feat.values' %in% which.plots){
      l_plots$feat.values <- p_feat_values
   }

   cowplot::plot_grid(
      plotlist=l_plots, nrow=1, align='h', axis='b',
      rel_widths=rel.widths
   )
}
