## Init ================================
options(stringsAsFactors=F)

WRITE_OUTPUT <- FALSE

## Paths --------------------------------
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')

## Dependencies --------------------------------
#devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
#library(matrixStats)

## Load data ================================
## Model --------------------------------
model_dir <- paste0(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/'))
out_dir <- paste0(model_dir,'/report/')
dir.create(out_dir, showWarnings=F)

metadata <- read.delim(paste0(model_dir,'/features/metadata_training.txt'))

features <- cacheAndReadData(paste0(model_dir,'/features/features_allSamples.rds'), overwrite=T)

anonymizeSampleId <- function(sample.id){
   #sample.id=pred_summ$sample
   #cohort <- metadata$cohort[match(sample.id, metadata$sample_id)]
   #is_hmf_sample <- cohort=='HMF'
   is_hmf_sample <- grepl('^(CPCT|WIDE|DRUP)',sample.id)

   out <- sample.id
   out[is_hmf_sample] <- metadata$sample_id_2[match(sample.id[is_hmf_sample], metadata$sample_id)]

   return(out)
}
# if(F){
#    features_exp <- features
#    rownames(features_exp) <- anonymizeSampleId(rownames(features_exp))
#    write.table(
#       features_exp, gzfile('HMF_PCAWG.cuplr_features.txt.gz', compression=9),
#       sep='\t',quote=F
#    )
# }

model <- cacheAndReadData(paste0(model_dir,'/final/model.rds'), overwrite=T)


## Pred reports --------------------------------
## CV set
gatherCvOutput(
   cv.out.dir=paste0(model_dir,'/cv_out/'),
   out.path=paste0(model_dir,'/test_set_report.rds')
)
pred_reports.cv <- cacheAndReadData(paste0(model_dir,'/test_set_report.rds'), overwrite=T)

## Val set samples
val_samples <- list(
   holdout=subset(metadata, is_holdout_sample, sample_id, drop=T),
   CUP=subset(metadata, cancer_type=='Unknown' & !is_blacklisted_sample, sample_id, drop=T),
   misc_ct=subset(metadata, !is_holdout_sample & !is_training_sample & !is_blacklisted_sample & cancer_type!='Unknown', sample_id, drop=T)
)

val_samples$excluded <- (function(){
   already_incl_samples <- c(
      unlist(val_samples, use.names=F),
      rownames(pred_reports.cv$prob)
   )
   sel_samples <- metadata$sample_id[ !(metadata$sample_id %in% already_incl_samples) ]
   sel_samples[sel_samples %in% rownames(features)]
})()

## HMF and PCAWG
pred_reports <- list(CV=pred_reports.cv)
pred_reports <- c(
   pred_reports,
   lapply(val_samples, function(i){
      out <- predict(
         object=model,
         newdata=features[i,],
         calc.feat.contrib=T, type='report', top.n.pred.classes=NULL
      )
      out$class_actual <- as.factor(
         metadata$cancer_type[ match(rownames(out$prob), metadata$sample_id) ]
      )

      return(out)
   })
)

## Fix the class of the objects
pred_reports <- lapply(pred_reports, function(i){
   class(i) <- c('predReport',class(i))
   return(i)
})


## Prob calibration ================================
cal_curves <- probCalCurves(report=pred_reports$CV, method='isotonic', prob.prescale=T)
#probCal(pred_reports$CV$prob, cal_curves)

pred_reports <- lapply(pred_reports, function(i){
   i$prob_scaled <- probCal(i$prob, cal.curves=cal_curves)
   i$class_pred <- colnames(i$prob_scaled)[ max.col(i$prob_scaled) ]
   i$class_pred <- factor(i$class_pred, colnames(i$prob_scaled))
   return(i)
})

## Summary table --------------------------------
pred_summ <- do.call(rbind, lapply(names(pred_reports), function(i){
   df <- summary(pred_reports[[i]], top.n.classes=3, top.n.feat=5, top.n.feat.classes=3, prob.type='prob_scaled')
   data.frame(group=i, df)
}))
#View(subset(pred_summ, group=='CUP'))

insColAfter <- function(df, v, after, colname=NULL){
   if(is.character(after)){
      after <- which(colnames(df)==after)
   }

   df_l <- df[1:after]
   df_r <- df[(after+1):ncol(df)]

   df_new <- cbind(df_l, v, df_r)

   if(!is.null(colname)){
      colnames(df_new)[(after+1)] <- colname
   }

   return(df_new)
}

pred_summ <- insColAfter(
   pred_summ, v=anonymizeSampleId(pred_summ$sample), after='sample', colname='sample_anon'
)

# if(WRITE_OUTPUT){
#    write.table(pred_summ, paste0(out_dir,'/pred_summ.txt'), sep='\t', quote=F, row.names=F)
# }

## Plot --------------------------------
if(WRITE_OUTPUT){
   p_cal_curves <- probCalCurves(report=pred_reports$CV, method='isotonic', prob.prescale=T, output='plot', facet.ncol=5)
   pdf(paste0(out_dir,'/prob_calib_curves.pdf'), 11, 11)
   plot(p_cal_curves)
   dev.off()

   write.table(cal_curves, paste0(out_dir,'/prob_calib_curves.txt'), sep='\t', quote=F, row.names=F)
}

## Analyse prob calibration --------------------------------
if(WRITE_OUTPUT){

   df_reliability <- (function(){
      ##
      df1 <- with(
         pred_reports$CV,
         reliabilityStats(actual=class_actual, probs=prob, verbose=T)
      )
      #df1$class <- paste0(df1$class,', uncal.')

      ##
      df2 <- with(
         pred_reports$CV,
         reliabilityStats(actual=class_actual, probs=prob_scaled, verbose=T)
      )
      df2$class <- paste0(df2$class,' (cal.)')

      ##
      df <- rbind(df1, df2)

      df$class <- factor(
         df$class,
         sort(unique(as.character(df$class)))
      )
      return(df)
   })()

   p_reliability <- reliabilityPlot(stats=df_reliability, facet.ncol=6)
   p_reliability <- colorizeFacetStrips(p_reliability, colors=c('lightgrey','lightsteelblue'), return.gtable=F)

   pdf(paste0(out_dir,'/reliability.pdf'), 11, 13)
   plot(p_reliability)
   dev.off()
}

## Main plots ================================
## Confusion matrix --------------------------------
if(WRITE_OUTPUT){

   pdf(paste0(out_dir,'/confusion_heatmap.pdf'), 11, 10.5)

   plot(with(
      pred_reports$CV,
      confusionHeatmap2(
         actual=class_actual, probs=prob_scaled,
         rel.heights=c(100,11,27), plot.title='Cross-validation (scaled probs)'
      )
   ))

   plot(with(
      pred_reports$CV,
      confusionHeatmap2(
         actual=class_actual, probs=prob,
         rel.heights=c(100,11,27), plot.title='Cross-validation (raw probs)'
      )
   ))

   plot(with(
      pred_reports$holdout,
      confusionHeatmap2(
         actual=class_actual, probs=prob_scaled,
         rel.heights=c(100,11,27), plot.title='Holdout (scaled probs)'
      )
   ))

   plot(with(
      pred_reports$holdout,
      confusionHeatmap2(
         actual=class_actual, probs=prob,
         rel.heights=c(100,11,27), plot.title='Holdout (raw probs)'
      )
   ))

   dev.off()
}

## Feature importance --------------------------------
if(WRITE_OUTPUT){

   imp_exp <- t(model$imp)
   imp_exp <- data.frame(feature=rownames(imp_exp), imp_exp, row.names=NULL)

   ## --------------------------------
   topFeaturesWrapper <- function(m, plot.title=NULL, facet.ncol=5){
      p <-
         topFeatures(m, infer.feature.type=T, top.n=15, facet.ncol=facet.ncol, feature.type.colors=FEATURE_TYPE_COLORS) +
         guides(fill=guide_legend(
            direction='horizontal',label.position='top',title.position='bottom',title.hjust=0.5,nrow=1,
            label.theme=element_text(angle=90, hjust=0, vjust=0.5, size=9)
         )) +
         theme(
            legend.position='bottom',
            legend.spacing.x=unit(0,'pt')
         )

      if(!is.null(plot.title)){
         p <- p + ggtitle(plot.title)
      }

      return(p)
   }

   pdf(paste0(out_dir,'/imp_top.pdf'), 11, 13)
   plot( topFeaturesWrapper(model$imp,'Final model') )
   plot( topFeaturesWrapper(pred_reports$CV$imp,'CV') )
   dev.off()

   # ## --------------------------------
   # p_imp_max <- maxImpPerFeatureType(pred_reports$CV$imp) +
   #    ylab('Cancer type') +
   #    theme(
   #       plot.margin=unit(c(0.01,0.04,0.01,0.01),'npc'),
   #       axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
   #       axis.text.x.top=element_blank(),
   #       axis.ticks.x.top=element_blank()
   #    )
   #
   # p_imp_sel <- (function(){
   #    sel_cancer_types <- c(
   #       'Cervix','CNS_PiloAstro','HeadAndNeck_SG','Lung_NSC',
   #       'Prostate','Sarcoma_GIST','Sarcoma_Lipo','Skin_Other'
   #    )
   #    m_imp <- model$imp[sort(sel_cancer_types),]
   #
   #    topFeatures(m_imp, infer.feature.type=T, top.n=15, facet.ncol=2, feature.type.colors=FEATURE_TYPE_COLORS) +
   #       guides(fill=guide_legend(
   #          direction='horizontal',label.position='bottom',title.position='top',title.hjust=0.5,nrow=1,
   #          label.theme=element_text(angle=90, hjust=1, vjust=0.5, size=9)
   #       )) +
   #       theme(
   #          #plot.margin=unit(c(0.01,0.04,0.01,0.01),'npc'),
   #          legend.position='bottom',
   #          legend.spacing.x=unit(0,'pt')
   #       )
   # })()
   #
   # p_imp_combined <- cowplot::plot_grid(
   #    p_imp_max, p_imp_sel,
   #    align='h', axis='b', nrow=1, rel_widths=c(1, 0.8)
   # )
   #
   # pdf(paste0(out_dir,'/imp_max.pdf'), 11, 8.5)
   # plot(p_imp_combined)
   # dev.off()

   # ## --------------------------------
   # p_imp_sel_list <- (function(){
   #    sel_cancer_types <- c('Cervix','Sarcoma_Lipo','Skin_Melanoma')
   #    m_imp <- model$imp[sort(sel_cancer_types),]
   #    topFeatures(
   #       m_imp, infer.feature.type=T, top.n=10, as.list=T, drop.legend.levels=F,
   #       feature.type.colors=FEATURE_TYPE_COLORS)
   # })()
   #
   # pdf(paste0(out_dir,'/imp_list.pdf'), 4, 3)
   # for(i in p_imp_sel_list) { plot(i) }
   # dev.off()
}

## Write table for feature description --------------------------------
if(WRITE_OUTPUT){
   df <- data.frame(
      feature_type=NA,
      feature_name=colnames(model$imp)
   )
   df$feature_type <- sub('[.].+','',df$feature_name)

   # googlesheets4::sheet_write(
   #    df,
   #    'https://docs.google.com/spreadsheets/d/1e0gxxgy-jqwtM20Yqxp-UXcoS_OLjZqawWW44-e-rDw/edit?usp=sharing',
   #    sheet='raw'
   # )
}

## Difference between CV and holdout perf --------------------------------
if(WRITE_OUTPUT){

   p_diff_cv_holdout <- (function(){
      ## Calc accuracy
      df_cv <- with(pred_reports$CV, calcFracCorrect(class_actual, class_pred))
      df_cv$group <- 'CV'

      df_h <- with(pred_reports$holdout, calcFracCorrect(class_actual, class_pred))
      df_h$group <- 'holdout'

      df <- data.frame(
         class=df_cv$class,
         frac_correct.cv=df_cv$frac_correct,
         frac_correct.holdout=df_h$frac_correct,
         n_total.cv=df_cv$n_total,
         n_total.holdout=df_h$n_total
      )
      df$n_total <- df$n_total.cv + df$n_total.holdout
      df <- subset(df, class!='All')

      ## Calc diff between holdout and CV accuracy
      relDiff <- function(x,y){ abs(x-y) / pmax(x,y) }
      df$frac_correct.perc_diff <- relDiff(df$frac_correct.h, df$frac_correct.cv)

      ## Point labels
      df$label <- paste0(df$class,' (',df$n_total.holdout,')')

      ## Plot
      ggplot(df, aes(x=n_total, y=frac_correct.perc_diff)) +
         geom_point(color='red') +
         ggrepel::geom_text_repel(aes(label=label), size=2.5, segment.size=0.2, nudge_y=0.01, nudge_x=40) +
         ylab('Relative difference in accuracy\nbetween CV and holdout sets') +
         xlab('x-axis: Total no. of samples per cancer type\nLabel values: No. of samples per cancer type in holdout set') +
         theme_bw()
   })()

   pdf(paste0(out_dir,'/diff_cv_holdout.pdf'), 7, 6)
   plot(p_diff_cv_holdout)
   dev.off()
}

## Compare original vs scaled probs --------------------------------
if(WRITE_OUTPUT){
   p_prob_raw_vs_scaled <- (function(){

      ## Predictions as dataframe
      sel_pred_types <- c('CV','holdout')

      df <- rbind(
         do.call(rbind, lapply(sel_pred_types, function(i){
            out <- summary(pred_reports[[i]], top.n.classes=3, prob.type='prob')
            data.frame(prob_type='raw', pred_type=i, out)
         })),

         do.call(rbind, lapply(sel_pred_types, function(i){
            out <- summary(pred_reports[[i]], top.n.classes=3, prob.type='prob_scaled')
            data.frame(prob_type='scaled', pred_type=i, out)
         }))
      )
      df$pred_type <- factor(df$pred_type, sel_pred_types)

      ## Add pan-cancer group
      df <- rbind(
         within(df,{ actual_class <- 'All' }),
         df
      )
      df$actual_class <- as.factor(df$actual_class)

      ## Summary statistics
      agg <- aggregate(
         df$pred_correct==1,
         list(actual_class=df$actual_class, pred_type=df$pred_type, prob_type=df$prob_type),
         function(x){
            c(n_correct=sum(x), n_total=length(x))
         }
      )
      agg <- cbind(agg[-ncol(agg)], agg$x)
      agg$frac_correct <- agg$n_correct / agg$n_total

      ## Plot
      agg <- agg[order(agg$pred_type),]
      agg$group <- paste0(agg$pred_type,'.',agg$prob_type)
      agg$group <- factor(agg$group, unique(agg$group))

      ggplot(agg, aes(x=group, y=frac_correct)) +
         facet_wrap(~actual_class) +

         geom_bar(aes(fill=group), stat='identity', color='black', size=0.3) +
         scale_fill_brewer(palette='Paired', direction=-1) +

         geom_text(
            aes(label=paste0(round(frac_correct,2),' (',n_correct,'/',n_total,')')),
            y=0.03, angle=90, hjust=0, size=2.7
         ) +

         ylab('Fraction correctly predicted') +

         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
         )
   })()

   pdf(paste0(out_dir,'/prob_scaled_comparison.pdf'), 11, 8)
   plot(p_prob_raw_vs_scaled)
   dev.off()
}


## Compare with other classifiers --------------------------------
if(WRITE_OUTPUT){

   metric_names <- c('recall','precision','f1')

   perf_comparison <- (function(){

      ## Calc CUPLR performance --------------------------------
      calcPerfWrapper <- function(report){
         #report=pred_reports$CV

         require(mltoolkit)

         confusion <- confusionMatrix(
            predicted=oneHotEncode(report$class_pred),
            actual=report$class_actual,
            simplify=T
         )

         perf <- calcPerf(confusion, metrics=metric_names)
         perf$n_samples <- confusion$tp + confusion$fn

         return(perf)
      }

      perf_cuplr <- calcPerfWrapper(pred_reports$CV)

      ## Compare perfs --------------------------------
      ## Get performance comparison template
      perf <- openxlsx::read.xlsx(paste0(base_dir,'/CUPs_classifier/processed/other_classifiers/classifier_perfs.xlsx'))

      ## Add performance for CUPLR
      perf_split <- split(perf, perf$classifier_name)

      getCuplrPerf <- function(col.name){
         perf_cuplr[
            match(perf_split$CUPLR$cancer_type_raw, perf_cuplr$class),
            col.name
         ]
      }
      perf_split$CUPLR$n_samples <- getCuplrPerf('n_samples')

      for(i in metric_names){
         perf_split$CUPLR[[i]] <- getCuplrPerf(i)
      }

      ## Calculate performance across all samples
      microAvg <- function(metric.values, n.samples){
         weights <- n.samples/sum(n.samples)
         #weights <- weights[names(i)]
         weighted.mean(metric.values, weights)
      }

      perf_split <- lapply(perf_split, function(i){
         #i=perf_split[[1]]
         i_All <- i[1,]
         'All' -> i_All$cancer_type_raw -> i_All$cancer_type_subgroup -> i_All$cancer_type_group -> i_All$cancer_type_color_group
         i_All$n_samples <- sum(i$n_samples)
         i_All$recall <- microAvg(i$recall, i$n_samples)
         i_All$precision <- microAvg(i$precision, i$n_samples)
         i_All$f1 <- microAvg(i$f1, i$n_samples)
         rbind(i_All, i)
      })

      perf <- do.call(rbind, perf_split); rownames(perf) <- NULL

      ## Force character order
      perf <- data.frame(lapply(perf, function(i){
         if(is.numeric(i)){ return(i) }
         factor(i, unique(i))
      }))

      ## Add number for bar coloring
      perf_split <- split(perf, perf$cancer_type_group)
      perf_split <- lapply(perf_split, function(i){
         #i=perf_split[[1]]
         i$color_index <- as.integer(droplevels(i$cancer_type_color_group))
         return(i)
      })

      perf <- do.call(rbind, perf_split); rownames(perf) <- NULL
      perf$color_index <- as.factor(perf$color_index)

      ## Melt
      perf <- do.call(rbind, lapply(metric_names, function(i){
         #i='recall'
         cbind(
            perf[,!(colnames(perf) %in% metric_names)],
            metric=i,
            metric_value=perf[,i]
         )
      }))

      ## Value labels
      perf$label <- with(perf,{
         paste0(cancer_type_raw,' (',round(metric_value,2),')')
      })

      perf <- perf[order(perf$classifier_name),]
      return(perf)
   })()

   lp_perf_comparison <- lapply(metric_names, function(i){
      #i='recall'
      ggplot(subset(perf_comparison, metric==i), aes(x=classifier_name, y=metric_value, group=cancer_type_subgroup)) +
         #facet_grid(classifier_name~cancer_type_group, scales='free_x') +
         facet_wrap(cancer_type_group~., ncol=7) +

         ## Main bars
         geom_bar(
            aes(fill=color_index), stat='identity', position=position_dodge(width=0.9),
            color='black', size=0.3, show.legend=F
         ) +
         scale_fill_brewer(palette='Set2') +

         ## Raw cancer type label
         geom_text(
            aes(label=label, y=0.02), position=position_dodge(width=0.9),
            angle=90, hjust=0, vjust=0.5, size=2.7
         ) +

         ## Axis labels
         ylab(tools::toTitleCase(i)) +
         xlab('Classifier name') +

         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
         )
   })
   names(lp_perf_comparison) <- metric_names

   lp_perf_comparison$recall <- lp_perf_comparison$recall + ylab('Fraction of samples correctly predicted')

   pdf(paste0(out_dir,'/perf_comparison.pdf'), 12, 11)
   #plot(lp_perf_comparison[[1]])
   for(i in lp_perf_comparison){ plot(i) }
   dev.off()
}

## Prob of target class ================================
if(WRITE_OUTPUT){
   pdf(paste0(out_dir,'/prob_target_class.pdf'), 11, 8)

   with(
      pred_reports$CV,
      targetClassProb(
         actual=class_actual, probs=prob_scaled,
         output='plot.per_cutoff', cutoffs=seq(0,1,0.1), show.labels=T
      ) + ggtitle('CV (samples grouped by actual class)')
   )

   with(
      pred_reports$holdout,
      targetClassProb(
         actual=class_actual, probs=prob_scaled,
         output='plot.per_cutoff', cutoffs=seq(0,1,0.1), show.labels=T
      ) + ggtitle('Holdout (samples grouped by actual class)')
   )

   with(
      pred_reports$CUP,
      targetClassProb(
         actual=class_pred, probs=prob_scaled,
         output='plot.per_cutoff', cutoffs=seq(0,1,0.1), show.labels=T
      ) + ggtitle('CUPs (samples grouped by predicted class)')
   )

   dev.off()
}

## Plot top-n accuracy ================================
# if(WRITE_OUTPUT){
#    pdf(paste0(out_dir,'/topn_acc.pdf'), 11, 8)
#    plot( with(pred_reports$CV, topnAcc(actual=class_actual, probs=prob_scaled)) + ggtitle('CV') )
#    plot( with(pred_reports$holdout, topnAcc(actual=class_actual, probs=prob_scaled)) + ggtitle('holdout') )
#    dev.off()
# }

if(WRITE_OUTPUT){
   p_topn_acc <- (function(){
      df_cv <- with(pred_reports$CV, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))
      df_holdout <- with(pred_reports$holdout, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))

      df <- rbind(
         cbind(pred_type='CV', df_cv),
         cbind(pred_type='Holdout', df_holdout)
      )

      df$label <- with(df,paste0(
         round(frac,2),' (',correct,'/',total,')'
      ))

      ggplot(df, aes(x=pred_type, y=frac, group=class_num)) +
         facet_wrap(~actual) +
         geom_bar(
            aes(fill=class_num),
            stat='identity', position='dodge',
            color='black', size=0.3, width=0.85
         ) +
         scale_fill_gradient(
            low='#4390BC', high='#EFF6B9', guide='legend', breaks=unique(df$class_num),
            name='Top-N class'
         ) +

         geom_text(
            aes(label=label, y=0.03),
            position=position_dodge(width=0.85),
            size=2.7, angle=90, hjust=0
         ) +

         scale_y_continuous(name='Top-N accuracy', limits=c(0,1)) +

         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            panel.grid.major.x=element_blank(),
            axis.title.x=element_blank()
         )
   })()

   pdf(paste0(out_dir,'/topn_acc.pdf'), 11, 8)
   plot(p_topn_acc)
   dev.off()
}

## Patient report examples ================================
# View(subset(pred_summ, group=='holdout'))
if(WRITE_OUTPUT){

   ##
   patientReportWrapper <- function(report, sample.name){

      sample_name_anon <- anonymizeSampleId(sample.name)
      cancer_type <- metadata$cancer_type[match(sample.name, metadata$sample_id)]
      plot_title <- paste0(sample_name_anon,' (',cancer_type,')')

      with(report,{
         patientReport(
            probs=prob_scaled,
            feat.contrib=feat_contrib,
            sample.name=sample.name,
            plot.title=plot_title,
            prob.thres.min=1,
            rel.widths=c(1.3, 1, 1), drop.feature.type.levels=F
         )
      })
   }

   patient_reports <- list(
      patientReportWrapper(pred_reports$holdout,'DO35949'), ## CNS_PiloAstro
      patientReportWrapper(pred_reports$holdout, sample.name='XXXXXXXX') ## Lung_SC
   )
   patient_reports <- cowplot::plot_grid(plotlist=patient_reports, ncol=1, rel_heights=c(1,1.6))

   pdf(paste0(out_dir,'/patient_reports.pdf'), 10, 11)
   patient_reports
   dev.off()

}

## Feature avg per cancer type ================================
if(WRITE_OUTPUT){
   p_feat_avg_per_ct <- (function(){
      df <- model$feat_stats

      ## Select top features
      ranked_features <- colnames(model$imp)[ order(-matrixStats::colMaxs(model$imp)) ]
      df <- df[df$feature %in% ranked_features[1:200],]

      ## Fix rmd.Pancreas.1 order
      insValueBefore <- function(x, value, before){
         if(is.character(before)){
            before <- match(before, x)
         }

         x_left <- x[1:(before-1)]
         x_right <- x[before:length(x)]
         c(x_left, value, x_right)
      }
      feature_names <- unique(df$feature)
      feature_names <- insValueBefore(feature_names,'rmd.Pancreas.1','rmd.Pancreas_NET.1')
      feature_names <- feature_names[ -head(which(feature_names=='rmd.Pancreas.1'),1) ]

      feature_names <- insValueBefore(feature_names,'rmd.Sarcoma_Lipo.1','rmd.Sarcoma_Other.1')
      feature_names <- feature_names[ -tail(which(feature_names=='rmd.Sarcoma_Lipo.1'),1) ]

      df$feature <- factor(df$feature, feature_names)
      df <- df[order(df$feature),]
      df$feature <- factor(df$feature, unique(df$feature))

      ## Parse feature names
      df$feature_type <- sub('[.].+$','',df$feature)
      df$feature_type <- factor(df$feature_type, unique(df$feature_type))
      df$feature_subtype <- sub('^\\w+[.]','',df$feature)
      df$feature_subtype <- factor(df$feature_subtype, unique(df$feature_subtype))

      ## Log transform features with a wide range
      df$is_wide_range <- with(df,{
         max_all-min_all > 100 & avg_metric != 'prop'
      })

      df$avg_case.trans <- df$avg_case
      df$min_all.trans <- df$min_all
      df$max_all.trans <- df$max_all

      df <- within(df,{
         avg_case.trans[is_wide_range] <- log10(avg_case.trans[is_wide_range]+1)
         min_all.trans[is_wide_range] <- log10(min_all.trans[is_wide_range]+1)
         max_all.trans[is_wide_range] <- log10(max_all.trans[is_wide_range]+1)
      })

      ## Rescale features from 0 to 1 for heatmap plotting
      scale0to1 <- function(x, x.min, x.max){
         out <- (x-x.min)/(x.max-x.min)

         ## Clip out of bounds feature values to (0,1)
         out[out>1] <- 1
         out[out<0] <- 0

         return(out)
      }
      df$avg_case.scaled <- with(df, scale0to1(avg_case.trans, min_all.trans, max_all.trans))

      ## Plotting
      default_text_size <- 7
      plotFeatureAvg <- function(df){
         ggplot(df, aes(x=as.integer(class), y=feature_subtype)) +
            facet_grid(feature_type~., scales='free_y', space='free_y') +

            geom_tile(aes(fill=avg_case.scaled), color='grey70') +
            scale_fill_distiller(
               palette='Spectral', name='Min max scaled\navg. cohort value\n',
               guide=guide_colorbar(
                  frame.colour='black', ticks.colour='black', barheight=5, barwidth=1,
                  label.theme=element_text(size=default_text_size),
                  title.theme=element_text(size=default_text_size)
               )
            ) +
            #scale_color_gradient2() +
            scale_y_discrete(expand=c(0,0), limits=rev) +
            scale_x_continuous(
               expand=c(0,0), sec.axis=dup_axis(),
               breaks=1:length(levels(df$class)),
               labels=levels(df$class)
            ) +

            theme_bw() +
            theme(
               panel.grid=element_blank(),
               panel.spacing = unit(3,'pt'),
               strip.text.y=element_text(angle=0, size=default_text_size),
               #strip.text.y.left=element_text(angle=0, size=default_text_size),
               #strip.placement='outside',
               axis.text.x.bottom=element_text(angle=90, vjust=0.5, hjust=1, size=default_text_size),
               axis.text.x.top=element_text(angle=90, vjust=0.5, hjust=0, size=default_text_size),
               axis.text.y=element_text(size=default_text_size),
               axis.title=element_blank()
            )
      }

      left_feature_types <- c('rmd','sigs')
      plots <- list()
      plots$left <- plotFeatureAvg( subset(df, feature_type %in% left_feature_types) ) + theme(legend.position='none')
      plots$right <- plotFeatureAvg( subset(df, !(feature_type %in% left_feature_types)) )

      cowplot::plot_grid(plotlist=plots, rel_widths=c(1, 1.2))
   })()

   pdf(paste0(out_dir,'/feat_avg_per_ct.pdf'), 10.5, 10.5)
   plot(p_feat_avg_per_ct)
   dev.off()
}

## Impact of MSI, HRD ================================
if(WRITE_OUTPUT){

   devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/statsExtra/'))

   incorrect_pred_enr <- (function(){
      ##
      df <- subset(pred_summ, group %in% c('CV','holdout'))
      df$actual_class <- as.character(df$actual_class)
      df$is_correct_pred <- df$actual_class==df$pred_class.1

      ## Add MSI and HRD status
      df <- cbind(
         df,
         metadata[
            match(df$sample, metadata$sample_id),
            c('msi_status','hr_status','smoking_history','treatment_platinum','treatment_5FU')
         ]
      )

      ## Parse characters to bool
      df$has_msi <- df$msi_status=='MSI'

      df$has_hrd <- NA
      df$has_hrd[df$hr_status=='HR_deficient'] <- TRUE
      df$has_hrd[df$hr_status=='HR_proficient'] <- FALSE

      df$has_smoked <- NA
      df$has_smoked[df$smoking_history %in% c('Current','Former')] <- TRUE
      df$has_smoked[df$smoking_history=='Never'] <- FALSE

      df$treated_with_platinum <- as.logical(df$treatment_platinum)
      df$treated_with_5FU <- as.logical(df$treatment_5FU)

      ##
      incorrectPredEnr <- function(df, colname, var.name=colname, sort.by.pvalue=TRUE){
         #colname='treatment_5FU'

         tab <- table(
            actual_class=df$actual_class,
            is_correct_pred=df$is_correct_pred,
            var_true=df[[colname]]
         )

         ## Contingency tables
         conting <- cbind(
            tab[,,'TRUE'],
            tab[,,'FALSE']
         )
         colnames(conting) <- c('var_true.incorrect','var_true.correct','var_false.incorrect','var_false.correct')
         conting <- as.data.frame(conting)

         ## Initialize output
         out <- data.frame(
            class=rownames(conting),
            n_samples=as.integer(table(df$actual_class)),
            n_samples.var_not_na=rowSums(conting),
            conting,
            row.names=NULL
         )

         ## Proportions
         out$var_true.n_samples <- with(conting, var_true.incorrect+var_true.correct)
         out$var_false.n_samples <- with(conting, var_false.incorrect+var_false.correct)

         out$var_true.frac_incorrect <- with(out, var_true.incorrect/var_true.n_samples)
         out$var_false.frac_incorrect <- with(out, var_false.incorrect/var_false.n_samples)

         ## Fisher test
         out$fisher_pvalue <- fisherTest(conting, alternative='greater')

         if(!is.null(var.name)){
            out <- cbind(var_name=var.name, out)
         }

         if(sort.by.pvalue){
            out <- out[order(out$fisher_pvalue),]
         }

         return(out)
      }



      rbind(
         incorrectPredEnr(df, 'has_msi'),
         incorrectPredEnr(df, 'has_hrd'),

         incorrectPredEnr(df, 'has_smoked'),

         incorrectPredEnr(df, 'treated_with_platinum'),
         incorrectPredEnr(df, 'treated_with_5FU')

      )
   })()

   googlesheets4::sheet_write(
      incorrect_pred_enr,
      'https://docs.google.com/spreadsheets/d/1x18ONHLl8F3R3fEYCuJfxYAo6TN2QozEtWvgk26bks0/edit?usp=sharing',
      sheet='raw'
   )

   # (function(){
   #    df <- incorrect_pred_enr
   #
   #    df$class_string <- with(df, paste0(class,' (',n_samples_with_ddrd_status,'/',n_samples,')'))
   #
   #    df$value_string.var_true <- with(df, paste0(
   #       round(var_true.frac_incorrect, 2),' (',var_true.incorrect,'/',var_true.n_samples ,')'
   #    ))
   #
   #    df$value_string.var_false <- with(df, paste0(
   #       round(var_false.frac_incorrect, 2),' (',var_false.incorrect,'/',var_false.n_samples ,')'
   #    ))
   #
   #    fixed_cols <- c('ddrd_name','class_string')
   #    df1 <- df[,c(fixed_cols,'var_true.frac_incorrect','value_string.var_true')]; df1$is_ddrd <- TRUE
   #    df2 <- df[,c(fixed_cols,'var_false.frac_incorrect','value_string.var_false')]; df2$is_ddrd <- FALSE
   #    c('frac_incorrect','value_string') -> colnames(df1)[3:4] -> colnames(df2)[3:4]
   #
   #    df_melt <- rbind(df1, df2)
   #    df_melt$value_string[df_melt$frac_incorrect==0 || is.na(df_melt$frac_incorrect)] <- ''
   #
   #    ggplot(subset(df_melt, ddrd_name=='MSI'), aes(x=class_string, y=frac_incorrect, group=is_ddrd)) +
   #       coord_flip() +
   #       #facet_wrap(.~ddrd_name, scales='free_x') +
   #       geom_bar(aes(fill=is_ddrd), stat='identity', position='dodge', color='black', size=0.3) +
   #       geom_text(aes(label=value_string, y=0.001), position=position_dodge(width=1), size=2.5, hjust=0) +
   #       scale_x_discrete(limits=rev) +
   #       theme_bw() +
   #       theme()
   # })()
}

## Plot RMD profiles ================================
if(WRITE_OUTPUT){
   p_rmd_profiles <- plotRmdProfile(t(model$rmd_sig_profiles))
   pdf(paste0(out_dir,'/rmd_profiles.pdf'), 11, 3)
   for(i in p_rmd_profiles){ plot(i) }
   dev.off()
}


