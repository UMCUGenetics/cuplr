#' Make barplots from a matrix of feature importances (for multiclass classification)
#'
#' @param m A numeric matrix where columns represent the features and rows represent the classes
#' (i.e. the names of each binary random forest)
#' @param cv_out Alternative input to `m`. Cross validation output in the form of a list that has
#' 'imp' object (e.g. cv_out[[1]]$imp)
#' @param top.n Top number of features to show
#' @param feature.type.colors A character vector of color hex codes with names being the feature
#' types
#' @param infer.feature.type Determine the feature type based on the tag/prefix. Everything before
#' @param infer.feature.type.func Function used to infer feature type. Defaults to
#' `function(x){ gsub('[.].+$','',x) }`, which takes the string before the first dot as the tag
#' the first dot (.) is considered the feature tag
#' @param facet.ncol Number of facet columns. If 'auto', this value will be
#' `round(sqrt(n_classes))`
#' @param facet.rnow Number of facet rows
#' @param as.list If FALSE (default), a faceted ggplot is returned (one facet for each class). If
#' TRUE, each class is plotted as separate ggplot and a list of these is returned
#'
#' @return A ggplot object
#' @export
#'
plotTopFeatures <- function(
   m=NULL, cv_out=NULL, top.n=10,
   feature.type.colors=NULL,
   infer.feature.type=F, infer.feature.type.func=NULL,
   facet.ncol='auto', facet.nrow=NULL,
   as.list=F
){
   if(F){
      m=readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.11a.1_snv_contexts/imp.rds')
      infer.feature.type.func=function(x){ gsub('(^\\w\\[)|(\\]\\w$)','',x) }

      m=readRDS('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.11c.1_indel_contexts/imp.rds')
      infer.feature.type.func=NULL
   }

   ## Init --------------------------------
   require(ggplot2)

   if(!is.null(cv_out)){
      m <- aggregateMatrixList(lapply(cv_out,`[[`,'imp'), as.matrix=T)
   }

   df <- do.call(rbind, lapply(rownames(m), function(i){
      v <- sort(m[i,], decreasing=T)[1:top.n]
      data.frame(class=i, feature=names(v), value=v, index=1:top.n, row.names=NULL)
   }))
   df <- forceDfOrder(df)

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
   label_y_pos <- max(df$value) * 0.05

   ## Bar colors
   if(is.null(feature.type.colors)){
      # color_pal <- c(
      #    RColorBrewer::brewer.pal(12, 'Set3'),
      #    RColorBrewer::brewer.pal(9, 'Pastel1'),
      #    RColorBrewer::brewer.pal(8, 'Pastel2')
      # )

      color_pal <- c(
         "#8DD3C7","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F", ##Excl: "#D9D9D9","#FFFFB3"
         "#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2",
         "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC"
      )

      color_pal <- structure(
         color_pal[1:length(levels(df$feature_type))],
         names=levels(df$feature_type)
      )
   } else if(feature.type.colors=='auto'){
      color_pal <- c(
         sigs="#8DD3C7",
         kataegis="#FFFFB3",
         gene_amp="#FB8072",
         gene_def="#80B1D3",
         aneuploidy="#FDB462",
         #arm_loss="#FDB462",
         #arm_gain="#B3DE69",
         purple='#BC80BD',
         sv_types="#FCCDE5",
         viral_ins="#D9D9D9",
         fusion="#FFED6F",
         rep_elem="#CCEBC5",
         rmd="#BEBADA"
      )
   } else {
      color_pal <- feature.type.colors
   }

   ## No. of facet columns
   if(facet.ncol=='auto'){
      facet.ncol <- round(sqrt(nrow(m)))
   }

   ## Plot --------------------------------
   main <- function(pd){
      p <- ggplot(pd, aes(x=index, y=value))

      if(!as.list){
         p <- p + facet_wrap(~class, nrow=facet.nrow, ncol=facet.ncol)
      }

      p <- p +
         geom_bar(aes(fill=feature_type), stat='identity') +
         scale_fill_manual(values=color_pal, limits=names(color_pal)) +
         geom_text(aes(label=label, y=label_y_pos), angle=90, hjust=0, size=2.5) +
         labs(y='Feature importance', x='Index', fill='Feature type') +
         theme(
            panel.grid.minor=element_blank()
         )

      if(length(unique(pd$feature_type))==1 | !feature_tags_exist){
         p <- p + guides(fill=F)
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
#' Make heatmap from a matrix of feature importances (for multiclass classification)
#'
#' @param m A numeric matrix where columns represent the features and rows represent the labels
#' @param cv_out Alternative input to `m`. Cross validation output in the form of a list that has
#' 'imp' object (e.g. cv_out[[1]]$imp)
#' @param top.n Top number of features to show
#' @param min.imp Features below this importance value are excluded from the plot
#' @param invert.y Invert the y-axis?
#'  @param sort.features Order features alphabetically?
#'
#' @return A ggplot object
#' @export
#'
plotFeatureImpHeatmap <- function(
   m=NULL, cv_out=NULL, top.n=50, min.imp=NULL, invert.y=T, sort.features=F
){

   if(!is.null(cv_out)){
      m <- aggregateMatrixList(lapply(cv_out,`[[`,'imp'), as.matrix=T)
   }

   if(sort.features){
      m <- m[,order(colnames(m))]
   }

   df <- reshape2::melt(as.matrix(m))
   colnames(df) <- c('class','feature','value')

   if(!is.null(top.n)){
      df2 <- df[order(df$value, decreasing=T),]
      feature_whitelist <- unique(df2$feature)[1:top.n]
      df <- df[df$feature %in% feature_whitelist,]
   } else if(!is.null(min.imp)){
      feature_whitelist <- as.character(unique(df[df$value >= min.imp,'feature']))
      df <- df[df$feature %in% feature_whitelist,]
   }

   df <- forceDfOrder(df)
   df$index <- as.integer(df$feature)

   if(invert.y){
      df$class <- factor(df$class, rev(levels(df$class)))
   }

   n_sel_features <- length(unique(df$feature))
   xlabel <-
      if(ncol(m) < n_sel_features ){
         'Features'
      } else {
         paste0('Features (top ',n_sel_features,'/',ncol(m),')')
      }

   require(ggplot2)
   p <- ggplot(df, aes(y=class,x=index)) +
      geom_tile(aes(fill=value)) +
      scale_x_continuous(
         sec.axis=dup_axis(), breaks=unique(df$index), labels=levels(df$feature),
         expand=c(0,0)
      ) +
      scale_fill_distiller(palette='YlGnBu') +

      labs(y='Class',x=xlabel,fill='Feat.\nimp.') +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5)
      )

   if(!is.null(min.imp)){
      p <- p + xlab(sprintf('Features\n(min.imp>=%s in at least 1 class)', min.imp))
   }

   return(p)
}
