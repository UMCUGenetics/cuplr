#' Plot top-n accuracy
#'
#' @description For samples within each actual class, plot the top-n accuracy of said actual class.
#'
#' @param report A list with the objects with the names: prob, class_actual
#' @param actual A vector of the actual classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param top.n Number of top classes to show stats for
#' @param output If 'plot', returns a ggplot object, else 'values' returns a dataframe with the
#' values
#' @param show.all.classes.stats Show the panel with the stats for all samples (i.e. not split by
#' actual class)
#' @param max.prob Filter for samples <= than this probability
#' @param label.direction Can be 'horizontal' or 'vertical'. Direction of the labels indicating the
#' values for each bar.
#'
#' @return A ggplot object or a dataframe
#' @export
#'
topnAcc <- function(
   report=NULL, actual=NULL, probs=NULL, top.n=3, output=c('plot','values'),
   show.all.classes.stats=T, max.prob=NULL, label.direction='horizontal'
){
   # if(F){
   #    report <- pred_reports$CV
   #
   #    actual=report$class_actual
   #    probs=report$prob_scaled
   #
   #    top.n=3
   #    show.all.classes.stats=T
   #    max.prob=NULL
   #
   #    output='plot'
   # }

   ## Init --------------------------------
   require(ggplot2)

   if(!is.null(report)){
      actual <- report$class_actual
      probs <- report$prob
   }

   if(length(actual)!=nrow(probs)){
      stop('length(actual) and nrow(probs) must be equal')
   }

   if(!is.factor(actual)){
      warning('`actual` is not a factor. Factor levels will be set automatically')
      actual <- as.factor(actual)
   }

   output <- match.arg(output, c('plot','values'))
   label.direction <- match.arg(label.direction, c('vertical','horizontal'))

   ## --------------------------------
   uniq_classes <- colnames(probs)
   top_classes <- t(apply(probs,1,function(i){
      uniq_classes[order(i, decreasing=T)]
   }))

   ## Filter to top.n classes
   if(!is.null(top.n)){
      if(length(top.n)==1){
         top_n <- 1:top.n
      } else {
         top_n <- top.n
      }
   }

   top_correct <- sweep(top_classes[,top_n],1, actual, '==')

   topn_correct <- t(apply(top_correct, 1, function(i){
      which_true <- which(i)
      if(length(which_true)!=0){
         i[ which_true:length(i) ] <- TRUE
      }
      return(i)
   }))

   ## --------------------------------
   ## Convert data to long form
   df <- reshape2::melt(topn_correct)
   colnames(df) <- c('sample','class_num','is_topn_correct')
   #df$is_topn_correct <- reshape2::melt(topn_correct)$value

   ## Additional data
   sample_indexes <- match(df$sample, rownames(probs))
   df$actual <- actual[ sample_indexes ]
   df$max_prob <- matrixStats::rowMaxs(probs)[ sample_indexes ]
   if(!is.null(max.prob)){
      df <- df[df$max_prob<=max.prob,]
   }

   if(show.all.classes.stats){
      df <- rbind(
         within(df,{ actual <- 'All' }),
         df
      )
      df$actual <- factor(df$actual, c('All',uniq_classes))
   }

   ## Summary
   df_agg <- with(df,{
      aggregate(
         is_topn_correct,
         list(class_num=class_num, actual=actual),
         sum
      )
   })
   colnames(df_agg)[ncol(df_agg)] <- 'count'

   n_samples <- with(df, {
      out <- aggregate(sample, list(actual), function(x){ length(unique(x)) })
      structure(out$x, names=as.character(out[,1]))
   })

   df_agg$n_samples <- as.integer( n_samples[df_agg$actual] )
   df_agg$frac <- df_agg$count / df_agg$n_samples

   if(output=='values'){ return(df_agg) }

   ## Plot --------------------------------
   ## Labels
   df_agg$actual.n_samples <- with(df_agg,paste0(
      actual,' (',n_samples,')'
   ))

   if(label.direction=='horizontal'){
      label_angle <- 0
      label_hjust <- 0.5
      label_vjust <- 0
      df_agg$label <- with(df_agg, paste0(
         round(frac,2), '\n(',count,')'
      ))
   } else {
      label_angle <- 90
      label_hjust <- 0
      label_vjust <- 0.5
      df_agg$label <- with(df_agg, paste0(
         round(frac,2),' (',count,')'
      ))
   }

   ## Main
   ggplot(df_agg, aes(x=class_num, y=frac, group=1)) +
      facet_wrap(~actual.n_samples) +
      geom_bar(aes(fill=class_num), stat='identity', color='black', size=0.3, show.legend=F) +
      scale_fill_gradient(low='#4390BC', high='#EFF6B9') +

      geom_text(
         aes(label=label), y=0.05, size=2.7,
         angle=label_angle, hjust=label_hjust, vjust=label_vjust
      ) +

      scale_y_continuous(name='Top-n accuracy', limits=c(0,1)) +
      xlab('Top-n class') +

      theme_bw() +
      theme(
         panel.grid.minor=element_blank()
      )
}
