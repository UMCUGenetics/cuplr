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
#' @param top.n.class.probs The number of bars to show for the predicted class probabilities
#' @param top.n.class.features The number of panels to display that show the details for the top
#' predicted classes (including the exact probabilities, and feature contributions).
#' @param top.n.features The number features to show in the feature contribution plots in the
#' details panels
#' @param prob.thres.min The feature contribution panels will be shown for predicted classes with a
#' probability higher than this threshold
#' @param prob.thres.rel.diff The feature contribution panels will be shown for the top predicted
#' classes, as long as the relative difference between the 1st and 2nd classes, 2nd and 3rd classes,
#' etc, are higher than this threshold. Note that `prob.thres.min` overrides this threshold
#' @param gender.feature.name The name of the feature specifying the gender of each sample
#' @param prob.label.size A numeric value specifying the size of the probability labels
#' @param feature.label.size A numeric value specifying the size of the feature labels
#' @param level.of.detail An integer specifying which subplots to show. 1: Probabilities, 2:
#' Feature contributions, 3: Feature values. Specifying e.g. 3 will result in plots 1 and 2 also
#' being shown
#' @param rel.widths A numeric vector specifying the relative widths of the subplots
#'
#' @return A cowplot of ggplot objects
#' @export
#'
patientReport <- function(
   probs=NULL, feat.contrib=NULL,
   sample.name=NULL, plot.title=sample.name,
   top.n.class.probs=10,
   top.n.class.features=3, top.n.features=5,
   prob.thres.min=0.1, prob.thres.rel.diff=0.4,
   gender.feature.name='gender.gender',
   prob.label.size=3.5, feature.label.size=3.5,
   drop.feature.type.levels=TRUE,
   level.of.detail=3, rel.widths=c(1.3, 1, 1)
){

   if(F){
      ##
      sample.name='XXXXXXXX' ## 2 top classes
      sample.name='DO35949'

      probs=pred_reports$holdout$prob_scaled
      feat.contrib=pred_reports$holdout$feat_contrib
      plot.title=sample.name

      gender.feature.name='gender.gender';
      top.n.class.probs=10;
      top.n.class.features=3; top.n.features=5;
      prob.thres.min=0.1; prob.thres.rel.diff=0.4;
      prob.label.size=3.5; feature.label.size=3.5;
      drop.feature.type.levels=FALSE
      level.of.detail=3; rel.widths=c(1.3, 1, 1)
   }

   ## Init ================================
   require(ggplot2)
   require(ggrepel)

   ## LOD
   if(!(level.of.detail %in% 1:3)){
      stop('`level.of.detail` must be 1, 2 or 3')
   }

   ## Colors ================================
   ## Pred rank highlighting
   colors_dark  <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999") ## Set1
   colors_light <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2") ## Pastel1

   ## Probs ========================================================================================
   ## Checks  --------------------------------
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

   ## Select classes --------------------------------
   probs_top <- probs[order(-probs)][1:top.n.class.features]

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
      expr = {
         rle( lagged_rel_diff>=prob.thres.rel.diff )$lengths[1]
      },
      warning = function(w){ return(1) }
   )

   last_sel_class.min_prob <- tryCatch(
      expr = { max(which(probs_top>=prob.thres.min)) },
      warning = function(w){ return(1) }
   )

   last_sel_class <- max(last_sel_class.rel_diff, last_sel_class.min_prob)

   probs_top <- probs_top[ 1:last_sel_class ]

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
   top_n_class_probs <- min(ncol(probs), top.n.class.probs)

   p_probs <- ggplot(pd_probs[1:top_n_class_probs,], aes(x=class, y=prob)) +

      geom_segment(aes(x=class, xend=class, y=0, yend=prob, color=class), size=2, show.legend=F) +
      scale_color_manual(values=line_color) +

      geom_label(
         aes(fill=class, label=label, hjust=prob),
         size=prob.label.size, show.legend=F, label.padding=unit(0.15,'lines')
      ) +
      scale_fill_manual(values=label_fill) +

      scale_y_continuous(name='Cancer type probability', limits=c(0, 1), breaks=seq(0,1,0.2), expand=c(0.05, 0.05)) +
      scale_x_discrete(name='Cancer type', limits=rev) +

      coord_flip() +
      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_text(size=10),
         axis.text.x=element_text(size=10)
      )

   ## Add title
   if(!is.null(plot.title)){
      p_probs <- p_probs + ggtitle(plot.title)
   }

   ## Output
   if(level.of.detail==1){
      return(p_probs)
   }

   ## Feat contrib =================================================================================
   ## Checks --------------------------------
   if(!is.data.frame(feat.contrib)){
      stop('`feat.contrib` must be a long form dataframe')
   }
   feat.contrib <- subset(feat.contrib, sample==sample.name)

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

   ## Assign a unique feature tag for each class
   ## Ensures feature rank is preserved for each class
   pd_feat$feature_tag <- paste0(pd_feat$binary_rf, '::', pd_feat$feature)
   pd_feat$feature_tag <- factor(pd_feat$feature_tag, rev(unique(pd_feat$feature_tag)))

   ## Parse feature type
   pd_feat$feature_type <- sub('[.].*$','',pd_feat$feature)
   pd_feat$feature_type <- factor(pd_feat$feature_type, names(FEATURE_TYPE_COLORS))

   ## Make labels --------------------------------
   ## Indicate if feature value is higher/lower than other cancer type samples
   pd_feat$label_value <- ifelse(
      pd_feat$sample_value > pd_feat$avg_ctrl,
      'Higher',
      'Lower'
   )
   pd_feat$label_value[pd_feat$sample_value==pd_feat$avg_ctrl] <- 'Equal'

   ## Make labels
   pd_feat <- within(pd_feat,{
      ## Remove value for boolean features
      label_value[avg_metric=='prop'] <- ''

      ## Substitute gender boolean feature with male/female
      label_value[feature==gender.feature.name] <- 'Female'
      label_value[feature==gender.feature.name & sample_value==0] <- 'Male'
   })

   pd_feat$label <- with(pd_feat,{
      out <- paste0(' (',label_value,')')
      out[nchar(label_value)==0] <- ''
      out <- paste0(feature, out)
      return(out)
   })

   ## Plot --------------------------------
   div_pos <- (1:top.n.features)-0.5
   div_pos[1] <- 0

   p_contrib <- ggplot(pd_feat, aes(x=feature_tag, y=contrib)) +
      coord_flip() +
      facet_wrap(binary_rf~., scales='free_y', ncol=1, strip.position='left') +

      ## Separating lines
      geom_vline(
         mapping=aes(xintercept=stage(feature_tag, after_scale=div_pos)),
         size=0.2, linetype='dotted'
      ) +

      ## Main bars
      geom_bar(
         aes(fill=feature_type), stat='identity',
         color='black', size=0.2, width=0.8, show.legend=F
      ) +
      scale_fill_manual(values=FEATURE_TYPE_COLORS, name='Feature type', drop=drop.feature.type.levels) +

      ## Feature labels
      geom_text(
         aes(label=label), y=max(pd_feat$contrib)*0.01,
         size=feature.label.size, angle=0, hjust=0
      ) +
      ylab('Feature contribution\n\nIn brackets: Feature value\ncompared to other cancer types') +

      theme_bw() +
      theme(
         panel.grid.major.y=element_blank(),
         panel.grid.major.x=element_blank(),
         panel.grid.minor.x=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text=element_text(size=10, face='bold')
      )

   ## Edit strip colors
   strip_colors <- structure(
      colors_light[ 1:length(names(probs_top)) ],
      names=names(probs_top)
   )
   p_contrib <- colorizeFacetStrips(p_contrib, strip_colors)

   ## Output --------------------------------
   if(level.of.detail==2){
      return(cowplot::plot_grid(
         p_probs, p_contrib,
         nrow=1, align='h', axis='bt', rel_widths=rel.widths
      ))
   }

   ## Feat values ==================================================================================
   pd_featval <- pd_feat

   ## Calculate min/max feature value for each feature
   feat_value_types <- c('min_all','avg_ctrl','avg_case','sample_value','max_all')
   #feat_value_types <- c('avg_ctrl','sample_value')
   pd_featval$stat_min <- apply( pd_featval[,feat_value_types], 1, min )
   pd_featval$stat_max <- apply( pd_featval[,feat_value_types], 1 ,max )

   ## Convert to long form dataframe
   pd_featval <- do.call(rbind, lapply(feat_value_types, function(i){
      cbind(
         pd_featval[,c('binary_rf','feature','feature_tag','feature_rank','stat_min','stat_max','avg_metric')], ## Fixed columns
         value_type=i,
         value=pd_featval[,i] ## Variable (value) columns
      )
   }))

   ## Force order of plotting points
   pd_featval$value_type <- factor(pd_featval$value_type, feat_value_types)

   ## Log transform features with a wide range  --------------------------------
   pd_featval$is_wide_range <- with(pd_featval,{
      stat_max-stat_min > 100 & avg_metric != 'prop'
   })
   pd_featval$value.trans <- pd_featval$value
   pd_featval$stat_min.trans <- pd_featval$stat_min
   pd_featval$stat_max.trans <- pd_featval$stat_max

   pd_featval <- within(pd_featval,{
      value.trans[is_wide_range] <- log10(value.trans[is_wide_range]+1)
      stat_min.trans[is_wide_range] <- log10(stat_min.trans[is_wide_range]+1)
      stat_max.trans[is_wide_range] <- log10(stat_max.trans[is_wide_range]+1)
   })

   ## Scale feature values to arbitrary units ranging from 0 to 1  --------------------------------
   scale0to1 <- function(x, x.min, x.max){
      out <- (x-x.min)/(x.max-x.min)

      ## Clip out of bounds feature values to (0,1)
      out[out>1] <- 1
      out[out<0] <- 0

      return(out)
   }
   pd_featval$value.scaled <- with( pd_featval, scale0to1(value.trans, stat_min.trans, stat_max.trans) )

   #v <- featUnits(feature.name)
   #View(v[!grepl('^rmd',names(v))])

   ## Make pretty values --------------------------------
   pd_featval$value_label <- with(pd_featval,{
      #value=pd_featval$value
      out <- signif(value, 3)

      is_big_num <- value>100000
      out[is_big_num] <- formatC(value[is_big_num], format="e", digits=2)

      # is_small_num <- value>0 & value < 0.001
      # out[is_small_num] <- formatC(value[is_small_num], format="e", digits=2)

      return(out)
   })
   #0.385>0 & 0.385<0.001

   ## Subset dataframes for each plot component --------------------------------
   ## Dummy data frame to initialize ggplot
   pd_featval.init <- subset(
      pd_featval,
      value_type=='min_all',
      select=c(feature_tag, feature, binary_rf, value.scaled, is_wide_range)
   )

   featUnits <- function(feature.name){
      #feature.name=colnames(features)[-1]

      ## Occurrence: default
      feature_units <- rep('Occurrence',length(feature.name))

      feature_units[grepl('^rmd[.]',feature.name)] <- 'Prop. of SBSs'
      feature_units[grepl('^mut_load[.]',feature.name)] <- '# of mutations'
      feature_units[grepl('^chrom_arm[.]',feature.name)] <- 'CN diff vs genome CN'
      feature_units[feature.name==gender.feature.name] <- 'Prop. female'

      ##
      feature_units[grepl('^sigs[.]SBS',feature.name)] <- 'Prop. of SBSs'
      feature_units[grepl('^sigs[.]DBS',feature.name)] <- 'Prop. of DBSs'
      feature_units[grepl('^sigs[.]ID',feature.name)]  <- 'Prop. of indels'

      ##
      feature_units[
         grepl('^sv[.]LINEs$',feature.name) |
         grepl('^sv[.](DEL|DUP)_',feature.name)
      ] <- 'Prop. of SVs'

      feature_units[
         feature.name %in% c('sv.n_events','sv.complex.n_events','sv.double_minutes','sv.foldbacks')
      ] <- '# of events'

      feature_units[feature.name=='sv.complex.largest_cluster'] <- 'SVs in cluster'

      structure(feature_units, names=feature.name)
   }
   pd_featval.init$x_pos <- as.integer(pd_featval.init$feature_tag)

   pd_featval.init$value.units <- featUnits(pd_featval.init$feature)
   pd_featval.init <- within(pd_featval.init, {
      value.units[is_wide_range] <- paste0(value.units[is_wide_range],' (log scale)')
   })

   ## Segment between case and ctrl averages
   pd_featval.direction <- reshape2::dcast(
      data=subset(pd_featval, value_type %in% c('avg_ctrl','sample_value')),
      formula=binary_rf+feature_tag~value_type,
      value.var='value.scaled'
   )

   ## Range
   pd_featval.range <- subset(
      pd_featval,
      value_type %in% c('min_all','max_all'),
      c(feature_tag, binary_rf, value.scaled, value_label, value_type)
   )
   pd_featval.range$hjust <- 0
   pd_featval.range$hjust[pd_featval.range$value_type=='max_all'] <- 1

   ## Main values
   sel_value_types <- c('sample_value','avg_case','avg_ctrl')
   pd_featval.main <- subset(pd_featval, value_type %in% sel_value_types)
   pd_featval.main$value_type <- factor(pd_featval.main$value_type, sel_value_types)

   ## Plot --------------------------------
   color.avg_ctrl <- 'grey50'
   color.avg_case <- 'midnightblue'
   color.sample_value <- 'red2'

   p_featval <-
      ggplot(
         #subset(pd_featval, value_type %in% c('avg_ctrl','sample_value','min_all','max_all')),
         pd_featval.init,
         aes(x=feature_tag, y=value.scaled)
      ) +

      coord_flip() +
      facet_wrap(binary_rf~., scales='free_y', ncol=1, strip.position='left') +

      ## Separating lines
      geom_vline(
         mapping=aes(xintercept=stage(feature_tag, after_scale=div_pos)),
         size=0.2, linetype='dotted'
      ) +

      ## Min/max line
      geom_segment( ## segment
         mapping=aes(x=feature_tag, xend=feature_tag, y=0, yend=1), size=0.2
      ) +
      geom_text( ## Min/max labels
         data=pd_featval.range,
         mapping=aes(label=value_label, hjust=hjust), nudge_x=-0.2, size=2.7
      ) +
      geom_text( ## Unit label
         data=pd_featval.init,
         mapping=aes(label=value.units), y=0.5, nudge_x=-0.2, size=2.7
      ) +

      ## Main values
      geom_segment( ## Arrow between ctrl average and sample value
         data=pd_featval.direction,
         mapping=aes(x=feature_tag, xend=feature_tag, y=avg_ctrl, yend=sample_value),
         arrow=arrow(length=unit(0.15,'cm'), ends='last', type='closed'),
         color=color.sample_value, size=0.5
      ) +
      ggrepel::geom_text_repel( ## Values
         data=pd_featval.main,
         mapping=aes(label=value_label, color=value_type),
         direction='x', nudge_x=0.3,
         min.segment.length=0, segment.size=0.25, size=2.7, fontface='bold',
         show.legend=F
      ) +
      geom_point( ## Plot dummy points to force legend keys to be points
         data=pd_featval.main,
         mapping=aes(color=value_type),
         alpha=0
      ) +
      scale_color_manual(
         values=c(color.sample_value, color.avg_case, color.avg_ctrl),
         labels=c('Value in sample','Avg. in target cancer type','Avg. in other cancer types')
      ) +
      guides(color=guide_legend(override.aes=list(alpha=1, shape=15, size=5))) +

      ##
      scale_y_continuous(breaks=c(0,1), labels=c('Min','Max'), name='Feature value') +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         strip.background=element_blank(),
         strip.text=element_blank(),
         legend.title=element_blank(),
         legend.position='bottom',
         legend.direction='vertical',
         legend.margin=margin(0,0,0,0),
         legend.box.margin=margin(-3,-10,1,-10),
         legend.key.height=unit(10,'pt')
      )

   ## Output --------------------------------
   cowplot::plot_grid(
      p_probs, p_contrib, p_featval,
      nrow=1, align='h', axis='bt', rel_widths=rel.widths
   )
}

####################################################################################################
# ## Feat values ================================
# pd_featval <- pd_feat
#
# ## Calculate min/max feature value for each feature
# feat_value_types <- c('min_all','avg_ctrl','avg_case','sample_value','max_all')
# #feat_value_types <- c('avg_ctrl','sample_value')
# pd_featval$stat_min <- apply( pd_featval[,feat_value_types], 1, min )
# pd_featval$stat_max <- apply( pd_featval[,feat_value_types], 1 ,max )
#
# ## Convert to long form dataframe
# pd_featval <- do.call(rbind, lapply(feat_value_types, function(i){
#    cbind(
#       pd_featval[,c('binary_rf','feature','feature_label','feature_rank','stat_min','stat_max','avg_metric')], ## Fixed columns
#       value_type=i,
#       value=pd_featval[,i] ## Variable (value) columns
#    )
# }))
#
# ## Force order of plotting points
# pd_featval$value_type <- factor(pd_featval$value_type, feat_value_types)
#
# ## Log transform features with a wide range  --------------------------------
# pd_featval$is_wide_range <- with(pd_featval,{
#    stat_max-stat_min > 100 & avg_metric != 'prop'
# })
# pd_featval$value.trans <- pd_featval$value
# pd_featval$stat_min.trans <- pd_featval$stat_min
# pd_featval$stat_max.trans <- pd_featval$stat_max
#
# pd_featval <- within(pd_featval,{
#    value.trans[is_wide_range] <- log10(value.trans[is_wide_range]+1)
#    stat_min.trans[is_wide_range] <- log10(stat_min.trans[is_wide_range]+1)
#    stat_max.trans[is_wide_range] <- log10(stat_max.trans[is_wide_range]+1)
# })
#
# ## Scale feature values to arbitrary units ranging from 0 to 1  --------------------------------
# scale0to1 <- function(x, x.min, x.max){
#    out <- (x-x.min)/(x.max-x.min)
#
#    ## Clip out of bounds feature values to (0,1)
#    out[out>1] <- 1
#    out[out<0] <- 0
#
#    return(out)
# }
# pd_featval$value.scaled <- with( pd_featval, scale0to1(value.trans, stat_min.trans, stat_max.trans) )
#
# ## Add plotting values to data --------------------------------
# ## Make values human readable
# # formatLabels <- function(x, avg.metric, sig.digits=3){
# #    x <- signif(x, sig.digits)
# #    out <- sapply(x, function(i){
# #       if(i<1e3){ return(i) }
# #       if(i<1e6){ return(paste0(i/1e3,'K')) }
# #       if(i<1e9){ return(paste0(i/1e6,'M')) }
# #       if(i<1e12){ return(paste0(i/1e9,'T')) }
# #       return(i)
# #    })
# #
# #    out[avg.metric=='prop'] <- paste0( x[avg.metric=='prop']*100, '%' )
# #
# #    return(out)
# # }
#
# formatLabels <- function(x, avg.metric, sig.digits=2){
#    #x <- pd_featval$value
#    #out <- x
#    out <- formatC(pd_featval$value, format="e", digits=2)
#    out[x==0] <- 0
#    out[x>=1 & x<1e4] <- round(x[x>=1 & x<1e4],0)
#    out[x>1e-2 & x<1] <- signif(x[x>1e-2 & x<1], 2)
#
#    return(out)
# }
#
# pd_featval$value_label <- with(pd_featval, formatLabels(value, avg_metric, sig.digits=2))
#
# ## Dummy data frame for min/max segment --------------------------------
# pd_featval.min_max <- subset(pd_featval, value_type=='min_all', select=c(feature_label, binary_rf))
#
# ## Segment between case and ctrl averages --------------------------------
# pd_featval.direction <- reshape2::dcast(
#    data=subset(pd_featval, value_type %in% c('avg_ctrl','sample_value')),
#    formula=binary_rf+feature_label~value_type,
#    value.var='value.scaled'
# )
#
# ## Main values --------------------------------
# pd_featval <- subset(pd_featval, value_type %in% c('avg_ctrl','sample_value'))
#
# pd_featval$value_label.hjust <- as.integer(
#    subset(pd_featval, value_type=='sample_value', value, drop=T) <=
#       subset(pd_featval, value_type=='avg_ctrl', value, drop=T)
# )
#
# pd_featval <- within(
#    pd_featval,
#    value_label.hjust[value_type=='avg_ctrl'] <- 1 - value_label.hjust[value_type=='avg_ctrl']
# )
#
# ## Plot --------------------------------
# p_featval <- ggplot(pd_featval,aes(x=feature_label, y=value.scaled)) +
#
#    coord_flip() +
#    facet_wrap(binary_rf~., scales='free_y', ncol=1, strip.position='left') +
#
#    ## Min/max segment
#    geom_segment(
#       data=pd_featval.min_max,
#       mapping=aes(x=feature_label, xend=feature_label, y=0, yend=1),
#       color='lightgrey', size=0.5
#    ) +
#
#    ## Main points
#    geom_label(aes(label=value_label, color=value_type, hjust=value_label.hjust), size=3, show.legend=F) +
#    scale_color_manual(values=c(avg_ctrl='black',sample_value='indianred')) +
#
#    ## Segment between ctrl average and sample value
#    geom_segment(
#       data=pd_featval.direction,
#       mapping=aes(x=feature_label, xend=feature_label, y=avg_ctrl, yend=sample_value),
#       arrow = arrow(length=unit(0.2,"cm"), ends="last", type="closed"),
#       color='indianred'
#    ) +
#
#    scale_y_continuous(expand=c(0.1,0.1)) +
#
#    theme_bw() +
#    theme(
#       panel.grid=element_blank(),
#       axis.text.y=element_blank(),
#       axis.title.y=element_blank(),
#       axis.ticks.y=element_blank(),
#       strip.background=element_blank(),
#       strip.text=element_blank()
#    )
#
# cowplot::plot_grid(
#    p_probs, p_contrib, p_featval,
#    nrow=1, align='h', axis='bt', rel_widths=c(1.5, 1, 1)
# )
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ## Prep data --------------------------------
# feat_value_types <- data.frame(
#    matrix(c(
#       'min_all',      'black',     108, 3,   1,
#       'avg_ctrl',     'black',     108, 5,   1,
#       'avg_case',     'indianred', 16,  3,   1,
#       'sample_value', 'indianred', 21,  2.5, 0,
#       'max_all',      'black',     108, 3,   1
#    ),
#    dimnames=list(NULL, c('type','color','shape','size','alpha')),
#    ncol=5, byrow=T)
# )
#
# constructAesValues <- function(colname){
#    values <- feat_value_types[,colname]
#    if(grepl('^\\d',values[1])){
#       values <- as.numeric(values)
#    }
#    structure(
#       values,
#       names=feat_value_types$type
#    )
# }
#
# ## Calculate min/max feature value for each feature
# pd_feat$stat_min <- apply(pd_feat[,feat_value_types$type],1,min)
# pd_feat$stat_max <- apply(pd_feat[,feat_value_types$type],1,max)
#
# ## Conver to long form dataframe
# pd_feat_values <- lapply(feat_value_types$type, function(i){
#    cbind(
#       pd_feat[,c('binary_rf','feature','feature_label','feature_rank','stat_min','stat_max','avg_metric')], ## Fixed columns
#       value_type=i,
#       value=pd_feat[,i] ## Variable (value) columns
#    )
# })
# pd_feat_values <- do.call(rbind, pd_feat_values)
#
# ## Force order of plotting points
# pd_feat_values$value_type <- factor(pd_feat_values$value_type, feat_value_types$type)
#
# ## Log transform features with a wide range
# pd_feat_values$is_wide_range <- with(pd_feat_values,{
#    stat_max-stat_min > 100 & avg_metric != 'prop'
# })
# pd_feat_values$value.trans <- pd_feat_values$value
# pd_feat_values$stat_min.trans <- pd_feat_values$stat_min
# pd_feat_values$stat_max.trans <- pd_feat_values$stat_max
#
# pd_feat_values <- within(pd_feat_values,{
#    value.trans[is_wide_range] <- log10(value.trans[is_wide_range]+1)
#    stat_min.trans[is_wide_range] <- log10(stat_min.trans[is_wide_range]+1)
#    stat_max.trans[is_wide_range] <- log10(stat_max.trans[is_wide_range]+1)
# })
#
# ## Scale feature values to arbitrary units ranging from 0 to 1
# scale0to1 <- function(x, x.min, x.max){
#    out <- (x-x.min)/(x.max-x.min)
#
#    ## Clip out of bounds feature values to (0,1)
#    out[out>1] <- 1
#    out[out<0] <- 0
#
#    return(out)
# }
# pd_feat_values$value.scaled <- with( pd_feat_values, scale0to1(value.trans, stat_min.trans, stat_max.trans) )
#
# # ## Top to bottom feature order. 2 line spacing betwen features
# # pd_feat_values$feature_rank <- -pd_feat_values$feature_rank * 2
#
# pd_feat_values$feature_rank <- -pd_feat_values$feature_rank
#
# ## Add plotting values to data --------------------------------
# ## Make values human readable
# formatNumHuman <- function(x, avg.metric, sig.digits=3){
#    x <- signif(x, sig.digits)
#    out <- sapply(x, function(i){
#       if(i<1e3){ return(i) }
#       if(i<1e6){ return(paste0(i/1e3,'K')) }
#       if(i<1e9){ return(paste0(i/1e6,'M')) }
#       if(i<1e12){ return(paste0(i/1e9,'T')) }
#       return(i)
#    })
#
#    out[avg.metric=='prop'] <- paste0( x[avg.metric=='prop']*100, '%' )
#
#    return(out)
# }
# pd_feat_values$value_label <- with(pd_feat_values, formatNumHuman(value, avg_metric, sig.digits=2))
#
# ## Add hjust to ensure label doesnt go past min/max
# pd_feat_values$hjust <- 0.5
# pd_feat_values$hjust[pd_feat_values$value.scaled>=0.95] <- 1
# pd_feat_values$hjust[pd_feat_values$value.scaled<=0.05] <- 0
#
#
# # ggplot(pd_feat_values, aes(x=value_type, y=feature_label)) +
# #    facet_wrap(~binary_rf, scales='free_y', ncol=1) +
# #    geom_tile(aes(fill=value.scaled), color='black', alpha=0.8) +
# #    geom_text(aes(label=value_label)) +
# #    scale_fill_distiller(palette='Spectral') +
# #    theme_bw() +
# #    theme(
# #       panel.grid=element_blank()
# #    )
#
# ## Min/max value labels --------------------------------
# ## Get original min/max values
# pd_feat_values.range <- reshape2::dcast(
#    data=subset(pd_feat_values, value_type %in% c('min_all','max_all')),
#    formula=binary_rf+feature_rank+feature+is_wide_range+avg_metric~value_type,
#    value.var='value'
# )
# pd_feat_values.range <- reshape2::melt(
#    data=pd_feat_values.range,
#    measure.vars=c('min_all','max_all')
# )
#
# pd_feat_values.range$value_label <- formatNumHuman(
#    pd_feat_values.range$value,
#    pd_feat_values.range$avg_metric
# )
#
# ## Assign justification and offsets
# pd_feat_values.range$value.scaled <- -0.02
# pd_feat_values.range$hjust <- 0
# pd_feat_values.range <- within(pd_feat_values.range,{
#    value.scaled[variable=='max_all'] <- 1.02
#    hjust[variable=='min_all'] <- 1
# })
#
# ## Segment between case and ctrl averages --------------------------------
# pd_feat_values.direction <- reshape2::dcast(
#    data=subset(pd_feat_values, value_type %in% c('avg_case','avg_ctrl')),
#    formula=binary_rf+feature_rank~value_type,
#    value.var='value.scaled'
# )
# pd_feat_values.direction <- within(pd_feat_values.direction,{
#    min <- pmin(avg_case, avg_ctrl)
#    max <- pmax(avg_case, avg_ctrl)
# })
#
# ## Plot --------------------------------
# legend_labels <- c('Avg. in other cancer types','Avg. in target cancer type','Value in patient')
#
# p_feat_values <- ggplot(pd_feat_values, aes(x=feature_rank, y=value.scaled)) +
#    coord_flip() +
#    facet_wrap(binary_rf~., scales='free_y', ncol=1, strip.position='left') +
#
#    ## Min/max segment
#    geom_segment(
#       data=pd_feat_values.range,
#       mapping=aes(x=feature_rank, xend=feature_rank, y=0, yend=1),
#       color='grey', size=0.2
#    ) +
#
#    # ## Min/max labels
#    # geom_text(
#    #    data=pd_feat_values.range,
#    #    mapping=aes(x=feature_rank, y=value.scaled, label=value_label, hjust=hjust),
#    #    vjust=0.4, size=3.5, color='grey60', fontface='bold'
#    # ) +
#
#    # ## Feature label
#    # geom_text(
#    #    data=pd_feat_values.range,
# #    mapping=aes(x=feature_rank+0.8, y=0, label=feature),
# #    hjust=0, size=3, color='grey35'
# # ) +
#
# ## Segment between case and ctrl averages
# geom_segment(
#    data=pd_feat_values.direction,
#    mapping=aes(x=feature_rank, xend=feature_rank, y=min, yend=max),
#    color='indianred'
# ) +
#
#    ## Feature stats points. Show value for patient.
#    geom_point(
#       data=subset(pd_feat_values, value_type %in% c('avg_ctrl','avg_case','sample_value')),
#       aes(color=value_type, shape=value_type, size=value_type, alpha=value_type),
#       fill='white'
#    ) +
#    geom_label(
#       data=subset(pd_feat_values, value_type %in% c('sample_value')),
#       aes(color=value_type, label=value_label, hjust=hjust),
#       size=3.5, label.padding=unit(0.15,'lines'), show.legend=F
#    ) +
#    scale_color_manual(values=constructAesValues('color'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
#    scale_shape_manual(values=constructAesValues('shape'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
#    scale_size_manual (values=constructAesValues('size'), labels=legend_labels, guide=guide_legend(reverse=TRUE)) +
#    scale_alpha_manual(values=constructAesValues('alpha'), labels=legend_labels, guide=guide_legend(reverse=TRUE, override.aes=list(alpha=NA))) + ## alpha=NA: so that point on legend is drawn but point on plot not drawn
#
#    ##
#    scale_x_discrete() +
#    scale_y_continuous(
#       name='Feature value (variable units)',
#       breaks=c(0, 1),
#       labels=c('Min','Max')
#    ) +
#
#    theme_bw() +
#    theme(
#       #panel.grid.major.x=element_blank(),
#       panel.grid.minor.x=element_blank(),
#       axis.title.y=element_blank(),
#       axis.text.y=element_blank(),
#       axis.ticks.y=element_blank(),
#       legend.position='bottom',
#       legend.direction='vertical',
#       legend.margin=margin(0,0,0,0),
#       legend.box.margin=margin(-10,-10, 0,-10),
#       legend.key.height=unit(10,'pt'),
#       legend.title=element_blank(),
#       strip.background=element_blank(),
#       strip.text=element_blank()
#    )
#
# # cowplot::plot_grid(
# #    p_probs, p_contrib, p_feat_values,
# #    nrow=1, align='h', axis='bt', rel_widths=c(1.5, 1, 0.9)
# # )
#
# ## Add blank title to make sure plots are aligned
# if(!is.null(plot.title)){
#    p_feat_values <- p_feat_values + ggtitle(' ') ## Add blank title to make sure plots are aligned
# }
#
# ## Add strip colors
# p_feat_values <- colorizeFacetStrips(p_feat_values, strip_colors)
#
#
# ## Combine plots ================================
# #which.plots=c('probs','feat.contribs','feat.values')
#
# l_plots <- list()
# if('probs' %in% which.plots){
#    l_plots$probs <- p_probs
# }
# if('feat.contribs' %in% which.plots){
#    l_plots$feat.contribs <- p_contrib
# }
# if('feat.values' %in% which.plots){
#    l_plots$feat.values <- p_feat_values
# }
# l_plots <- l_plots[which.plots]
#
# cowplot::plot_grid(
#    plotlist=l_plots, nrow=1, align='h', axis='b',
#    rel_widths=rel.widths
# )
