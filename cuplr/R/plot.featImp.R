FEATURE_TYPE_COLORS <- c(
   rmd="#BEBADA", ## Set3:3 Lavender
   sigs="#8DD3C7", ## Set3:1 Teal
   mut_load="#DAB986", ## Misc Brown
   gene="#FCCDE5", ## Set3:8 Pink
   chrom_arm="#80B1D3", ## Set3:5 Blue
   genome="lightblue1",
   gender="#FB8072", ## Set3:4 Red
   sv="#FFFFB3", ## Set3:2 Pale yellow
   fusion="#FDB462", ## Set3:6 Light orange
   viral_ins="#B3DE69" ## Set3:7 Lime
)

####################################################################################################
#' Plot feature importance
#'
#' @rdname plot.featImp
#'
#' @description topFeatures(): Plot sorted feature importance barplots per class from a matrix of
#' feature importances. maxImpPerFeatureType(): Plot the max importance per feature type as a
#' heatmap
#'
#' @param m A numeric matrix where columns represent the features and rows represent the classes
#' (i.e. the names of each binary random forest)
#' @param top.n Top number of features to show
#' @param feature.type.colors A character vector of color hex codes with names being the feature
#' types
#' @param infer.feature.type Determine the feature type based on the tag/prefix. Everything before
#' @param infer.feature.type.func Function used to infer feature type. Defaults to
#' `function(x){ gsub('[.].+$','',x) }`, which takes the string before the first dot as the tag
#' @param facet.ncol Number of facet columns. If 'auto', this value will be
#' `round(sqrt(n_classes))`
#' @param facet.nrow Number of facet rows
#' @param axis.x.breaks A numeric vector indicating the x-axis breaks
#' @param hide.legend Hide legend?
#' @param drop.legend.levels Remove legend keys for feature types that don't exist in the plot?
#' @param as.list If FALSE (default), a faceted ggplot is returned (one facet for each class). If
#' TRUE, each class is plotted as separate ggplot and a list of these is returned
#'
#' @return A ggplot object
#' @export
#'
topFeatures <- function(
   m=NULL, top.n=10,
   feature.type.colors=NULL,
   infer.feature.type=F, infer.feature.type.func=NULL,
   facet.ncol=NULL, facet.nrow=NULL,
   axis.x.breaks=NULL, hide.legend=F, drop.legend.levels=T,
   as.list=F
){
   if(F){
      m=pred_reports$CV$imp
      top.n=10
      feature.type.colors=NULL
      infer.feature.type=T
      infer.feature.type.func=NULL
      facet.ncol='auto'
      facet.nrow=NULL
      as.list=F
   }

   ## Init --------------------------------
   require(ggplot2)

   df <- do.call(rbind, lapply(rownames(m), function(i){
      v <- sort(m[i,], decreasing=T)[1:top.n]
      data.frame(class=i, feature=names(v), value=v, index=1:top.n, row.names=NULL)
   }))
   df <- as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))

   ## Infer feature types --------------------------------
   feature_tags_exist <- all(grepl('^\\w+[.]',colnames(m)))
   if(infer.feature.type=='force' | (infer.feature.type==T & feature_tags_exist)){

      if(is.null(infer.feature.type.func)){
         infer.feature.type.func <- function(x){ gsub('[.].+$','',x) }
      }

      df$feature_type <- factor(
         infer.feature.type.func(df$feature),
         levels=unique(infer.feature.type.func(colnames(m)))
      )
   } else {
      df$feature_type <- 'none'
   }

   ## Plot params --------------------------------
   ## Bar labels
   df$label <- as.character(df$feature)
   df$label[df$value<=0] <- ''

   label_ypos <- with(df,{
      class_maxes <- aggregate(value, list(class), max)
      class_maxes <- structure(class_maxes$x, names=as.character(class_maxes$Group.1))
      class_maxes * 0.05
   })
   df$label_ypos <- label_ypos[ as.character(df$class) ]

   ## Bar colors
   if(is.null(feature.type.colors)){
      # color_pal <- c(
      #    RColorBrewer::brewer.pal(12, 'Set3'),
      #    RColorBrewer::brewer.pal(9, 'Pastel1'),
      #    RColorBrewer::brewer.pal(8, 'Pastel2')
      # )

      color_pal <- c(
         ## Muted colors
         "#8DD3C7", ## Set3:1 Teal
         "#BEBADA", ## Set3:3 Lavender
         "#DAB986", ## Misc Brown
         "#D9D9D9", ## Set3:9 Grey
         "#FCCDE5", ## Set3:8 Pink

         ## Bright colors
         "#FDB462", ## Set3:6 Light orange
         "#80B1D3", ## Set3:5 Blue
         "#FB8072", ## Set3:4 Red
         "#B3DE69", ## Set3:7 Lime
         "#FFFFB3", ## Set3:2 Pale yellow

         ## Recycled colors, different shade
         "#CCEBC5", ## Set3:11 Soap green
         "#BC80BD", ## Set3:10 Violet
         "#E0C49A",  ## Set2:7  Light brown
         "#F9DAAD", ## Pastel1:5 Pale orange
         "#FFED6F" ## Set3:12 Yellow
      )

      feature_type_counts <- sort(table(df$feature_type), decreasing=T)
      color_pal <- structure(
         color_pal[1:length(feature_type_counts)],
         names=names(feature_type_counts)
      )

   } else {
      color_pal <- feature.type.colors
   }

   ## Remove non-existing feature types
   if(drop.legend.levels){
      color_pal <- color_pal[names(color_pal) %in% df$feature_type]
   }

   ## Plot --------------------------------
   main <- function(pd){
      #pd=df
      p <- ggplot(pd, aes(x=index, y=value))

      if(!as.list){
         p <- p + facet_wrap(~class, nrow=facet.nrow, ncol=facet.ncol)
      }

      p <- p +
         geom_bar(aes(fill=feature_type), stat='identity', width=1, size=0.25, color='grey25') +
         scale_fill_manual(values=color_pal, limits=names(color_pal), drop=drop.legend.levels) +
         geom_text(aes(label=label, y=label_ypos), angle=90, hjust=0, size=2.5, color='black') +
         labs(y='Feature importance (mean decrease in accuracy)', x='Rank', fill='Feature type') +
         theme_bw() +
         theme(
            panel.grid.minor=element_blank()
         )

      if(!is.null(axis.x.breaks)){
         axis.x.breaks
      } else if(top.n<=10){
         axis_x_breaks <- seq(0,top.n,2)
      } else {
         axis_x_breaks <- waiver()
      }
      p <- p + scale_x_continuous(breaks=axis_x_breaks)

      if(length(unique(pd$feature_type))==1 | !feature_tags_exist){
         p <- p + guides(fill=F)
      }

      if(hide.legend){
         p <- p + theme(legend.position='none')
      }

      return(p)
   }

   if(!as.list){
      out <-  main(df)
   } else {
      df_split <- split(df, df$class)
      out <- lapply(names(df_split), function(i){
         main(df_split[[i]]) + ggtitle(i)
      })
   }

   return(out)
}

####################################################################################################
#' @rdname plot.featImp
#' @export
#'
maxImpPerFeatureType <- function(m, infer.feature.type.func=NULL){
   if(F){
      m=pred_reports$CV$imp
   }

   ## Init --------------------------------
   m <- as.matrix(m)
   if(!is.numeric(m)){ stop('`m` must be a numeric matrix or dataframe') }

   ## Infer feature types --------------------------------
   feature_tags_exist <- all(grepl('^\\w+[.]',colnames(m)))
   if(!feature_tags_exist){ stop('Feature names must be in the form: {feature tag}.{feature name}') }

   if(is.null(infer.feature.type.func)){
      infer.feature.type.func <- function(x){ gsub('[.].+$','',x) }
   }

   feature_types <- infer.feature.type.func(colnames(m))
   levels(feature_types) <- unique(feature_types)

   ## Stats --------------------------------
   max_imps <- lapply(
      split(as.data.frame(t(m)),feature_types),
      function(i){ matrixStats::colMaxs(as.matrix(i)) }
   )
   max_imps <- do.call(cbind, max_imps)
   rownames(max_imps) <- rownames(m)
   max_imps <- max_imps[,order(matrixStats::colMedians(max_imps), decreasing=T)]

   ## Plot --------------------------------
   df <- reshape2::melt(max_imps)
   colnames(df) <- c('class','feature_type','max_imp')
   df$feature_type <- factor(df$feature_type, unique(df$feature_type))
   df$class <- factor(df$class, unique(df$class))
   df$x_index <- as.integer(df$feature_type)

   ggplot(df, aes(x=x_index, y=class)) +
      geom_tile(aes(fill=max_imp), color='black', size=0.3) +
      geom_text(aes(label=round(max_imp,2)), size=2.7) +
      scale_fill_distiller(
         name='Max importance per feature type\n(mean decrease in accuracy)', palette='Spectral',
         guide=guide_colorbar(frame.colour='black', ticks.colour='black', title.vjust=1, title.position='right')
      ) +

      scale_x_continuous(
         name='Feature type', expand=c(0,0),
         breaks=unique(df$x_index), labels=levels(df$feature_type), sec.axis=dup_axis()
      ) +
      scale_y_discrete(name='Class', expand=c(0,0), limits=rev) +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.title.x.top=element_blank(),
         axis.text.x.top=element_text(angle=45, hjust=0),
         axis.text.x.bottom=element_text(angle=45, hjust=1),
         legend.position='bottom'
      )
}
