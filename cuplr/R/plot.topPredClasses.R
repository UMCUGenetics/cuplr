#' Plot top predicted classes
#'
#' @param actual A factor vector of the actual classes
#' @param predicted A factor vector of the predicted classes
#' @param top.n The top predicted classes to show
#' @param output If 'plot', returns a ggplot object, else 'raw' returns the raw data
#' @param show.all.classes.stats Show stats for all samples (i.e. not split by actual class)
#' @param sort.classes Sort classes by the predicted class?
#' @param show.labels Show labels indicating the numerical values at each point?
#'
#' @return A ggplot object
#' @export
#'
topPredClasses <- function(
   actual, predicted=NULL, probs=NULL, top.n=2, output=c('plot','raw'),
   show.all.classes.stats=T, sort.classes=F, show.labels=F
){

   if(F){
      actual=pred_reports$CV$class_actual
      predicted=pred_reports$CV$class_pred
   }

   ## Init --------------------------------
   require(ggplot2)

   if(!is.null(probs)){
      predicted <- factor(
         colnames(probs)[ max.col(probs) ],
         colnames(probs)
      )
   }

   if(length(actual)!=length(predicted)){
      stop('length(actual) and nrow(probs) must be equal')
   }

   output <- match.arg(output, c('plot','raw'))

   ## Format factor levels --------------------------------
   pred_levels <- levels(predicted)

   uniq_responses <- sort(unique(
      c(
         as.character(actual),
         as.character(predicted)
      )
   ))

   actual <- factor(actual, uniq_responses)
   predicted <- factor(predicted, uniq_responses)

   ## Calc stats --------------------------------
   ##
   tab <- table(predicted, actual)
   n_samples <- table(actual)
   n_samples <- structure(as.integer(n_samples), names=names(n_samples))

   df <- do.call(rbind, lapply(colnames(tab), function(i){
      #i='Bile_Gallbladder'
      v <- sort(tab[,i], decreasing=T)
      data.frame(
         actual=i,
         predicted=names(v),
         pred_rank=1:length(v),
         n_pred=v,
         n_samples=n_samples[i],
         prop_pred=v/n_samples[i],
         row.names=NULL
      )
   }))

   ## Filter top pred classes
   if(!is.null(top.n)){
      df <- df[df$pred_rank <= top.n,]
   }

   if(show.all.classes.stats){
      df <- rbind(
         (function(){
            df1 <- df[1,]
            df1$actual <- 'All'
            df1$predicted <- 'All'
            df1$n_pred <- sum(diag(tab))
            df1$n_samples <- length(actual)
            df1$prop_pred <- df1$n_pred / df1$n_samples
            return(df1)
         })(),
         df
      )
      df$actual <- factor(df$actual, unique(df$actual))
   }

   ## Sort classes by accuracy
   class_order <- levels(actual)
   if(sort.classes){
      df_top1 <- df[df$pred_rank==1,]
      df_top1 <- df_top1[order(df_top1$prop_pred, decreasing=T),]
      class_order <- as.character(df_top1$actual)
   }

   if(show.all.classes.stats){
      class_order <- unique(c('All',class_order))
   }

   df$actual <- factor(df$actual, class_order)

   df <- df[order(df$actual),]

   if(output=='raw'){ return(df) }

   ## Plot --------------------------------
   ## Labels
   df$label <- with(df, paste0(
      predicted,': ',
      round(prop_pred,2),' (',n_pred,')'
   ))

   df$actual.n_samples <- with(df, paste0(
      actual,' (',n_samples,')'
   ))
   df$actual.n_samples <- factor(df$actual.n_samples, unique(df$actual.n_samples))

   df$pred_rank_2 <- paste0('Pred. class #',df$pred_rank)
   df$pred_rank_2 <- factor(df$pred_rank_2, unique(df$pred_rank_2))

   ## Main
   p <- ggplot(df, aes(x=actual.n_samples, y=prop_pred))

   if(top.n>1){
      p <- facet_grid(pred_rank_2~.)
   }

   p <- p +
      geom_bar(aes(fill=pred_rank_2, color=pred_rank_2), stat='identity', size=0.25, show.legend=F) +
      scale_fill_brewer(palette='Pastel1') +
      scale_color_brewer(palette='Set1') +
      scale_x_discrete(position='top', name='Actual class; total # of samples') +
      labs(y='Proportion of samples') +

      theme_bw() +
      theme(
         panel.grid.minor.x=element_blank(),
         panel.grid.minor.y=element_blank(),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         strip.text.y=element_text(angle=90)
      )

   if(show.labels){
      p <- p + geom_text(aes(label=label), y=0.01, angle=90, hjust=0, size=2.7)
   }

   return(p)
}
