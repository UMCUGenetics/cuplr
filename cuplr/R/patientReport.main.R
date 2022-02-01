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
   top.n.class.probs=NULL, top.n.class.features=3, top.n.features=5,
   prob.thres.min=0.1, prob.thres.rel.diff=0.4,
   gender.feature.name='gender.gender',
   drop.feature.type.levels=TRUE,
   level.of.detail=3, rel.widths=c(1.05, 1, 1)
){

   if(F){
      ##
      sample.name='XXXXXXXX' ## 2 top classes
      #sample.name='DO36039'

      probs=pred_reports$holdout$prob_scaled
      feat.contrib=pred_reports$holdout$feat_contrib
      plot.title=sample.name

      gender.feature.name='gender.gender';
      top.n.class.probs=NULL;
      top.n.class.features=3; top.n.features=5;
      prob.thres.min=0.1; prob.thres.rel.diff=0.4;
      prob.label.size=3.5; feature.label.size=3.5;
      drop.feature.type.levels=FALSE
      level.of.detail=3; rel.widths=c(1.05, 1, 1)
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

   pd_probs$prob_round <- round(pd_probs$prob, 2)
   pd_probs$label <- format(pd_probs$prob_round, nsmall=2)
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
   line_color[pd_probs$prob_round==0] <- 'grey50'

   label_fill <- structure(
      colors_light[color_indexes],
      names=as.character(pd_probs$class)
   )
   label_fill[pd_probs$prob_round==0] <- 'lightgrey'

   ## Main --------------------------------
   if(is.null(top.n.class.probs)){
      ## Show 10 probabilities per class with features shown
      ## More class features shown = more probabilities shown
      top_n_class_probs <- length(probs_top)*5
   } else {
      top_n_class_probs <- top.n.class.probs
   }
   top_n_class_probs <- min(length(probs), top_n_class_probs) ## Prevent `top_n_class_probs` from exceeding the total number of predicted classes

   ##
   pd_probs <- pd_probs[1:top_n_class_probs,]

   #row_spacing <- 2.5
   pd_probs$index <- nrow(pd_probs):1
   pd_probs$index <- pd_probs$index #* row_spacing
   text_offset <- 0.5
   bar_offset <- 0

   p_probs <-
      ggplot(pd_probs, aes(x=index, y=prob)) +

      ## Min/max line
      geom_segment(
         aes(x=index+bar_offset, xend=index+bar_offset, y=0, yend=1),
         size=0.5, color='grey'
      ) +

      ## Lollipop
      geom_segment(
         aes(x=index+bar_offset, xend=index+bar_offset, y=0, yend=prob, color=class),
         size=2.5, show.legend=F
      ) +
      geom_label(
         aes(fill=class, label=label, hjust=prob, x=index+bar_offset),
         size=3.5, show.legend=F, label.padding=unit(0.15,'lines')
      ) +

      ## Cancer type labels
      geom_text(aes(label=class, y=0, x=index+text_offset, color=class), hjust=0, vjust=1, show.legend=F) +

      ## Colors
      scale_color_manual(values=line_color) +
      scale_fill_manual(values=label_fill) +

      ## Axes
      scale_y_continuous(
         name='Cancer type probability',
         breaks=seq(0, 1, 0.2), limits=c(0, 1), expand=c(0.04, 0.04)
      ) +
      scale_x_continuous(
         name='Prediction rank',
         breaks=seq(1, max(pd_probs$index), 1),
         limits=c(min(pd_probs$index)-0.6, max(pd_probs$index)+0.4+text_offset),
         expand=c(0,0),
         labels=function(x){ max(x, na.rm=T)-x+1 }
      ) +
      coord_flip() +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
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
         angle=0, hjust=0
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
      feature_units[grepl('^chrom_arm[.]',feature.name)] <- 'Arm CN/genome CN ratio'
      feature_units[feature.name==gender.feature.name] <- 'Prop. female'
      feature_units[feature.name=='genome.diploid_proportion'] <- 'Prop. of genome'
      feature_units[feature.name=='genome.ploidy'] <- 'Genome ploidy'
      feature_units[feature.name=='sv.COMPLEX.largest_cluster'] <- '# of breakpoints in largest cluster'

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
      ggrepel::geom_text_repel( ## Values
         data=pd_featval.main,
         mapping=aes(label=value_label, color=value_type),
         direction='x', nudge_x=0.25,
         min.segment.length=0, segment.size=0.4, size=2.7, fontface='bold',
         show.legend=F
      ) +
      geom_point( ## Plot dummy points to force legend keys to be points
         data=pd_featval.main,
         mapping=aes(color=value_type),
         alpha=0
      ) +
      scale_color_manual(
         values=c(sample_value='red3', avg_case='lightcoral', avg_ctrl='dodgerblue4'),
         labels=c(sample_value='Value in sample', avg_case='Avg. in target cancer type', avg_ctrl='Avg. in all other samples')
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
