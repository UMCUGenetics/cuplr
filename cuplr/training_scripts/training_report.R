## Init ================================
options(stringsAsFactors=F)

WRITE_OUTPUT <- TRUE

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
model_dir <- paste0(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.18a_newLabels/'))
out_dir <- paste0(model_dir,'/report/')
dir.create(out_dir, showWarnings=F)

model <- cacheAndReadData(paste0(model_dir,'/final/model.rds'), overwrite=T)
features <- cacheAndReadData(paste0(model_dir,'/features/features_allSamples.rds'), overwrite=T)
metadata <- read.delim(paste0(model_dir,'/features/training_metadata.txt.gz'))
# with(
#    subset(metadata, in_train_set | in_val_set),
#    table(cohort)
# )

# xx <- predict(
#    object=model,
#    newdata=features[1:10,],
#    calc.feat.contrib=T, type='report', top.n.pred.classes=NULL
# )

## Pred reports --------------------------------
## CV set
gatherCvOutput(
   cv.out.dir=paste0(model_dir,'/cv_out/'),
   out.path=paste0(model_dir,'/test_set_report.rds')
)
pred_reports.cv <- cacheAndReadData(paste0(model_dir,'/test_set_report.rds'), overwrite=T)

## Val set samples
misc_ct <- metadata$cancer_type[
   !(metadata$cancer_type %in% names(model$ensemble)) &
   metadata$cancer_type!='Unknown'
]

val_samples <- list(
   holdout=subset(metadata, in_val_set, sample_id, drop=T),
   CUP=subset(metadata, cancer_type=='Unknown' & tumor_purity_rank==1, sample_id, drop=T),
   misc_ct=subset(metadata, cancer_type %in% misc_ct & tumor_purity_rank==1, sample_id, drop=T)
)

val_samples$excluded <- metadata$sample_id[
   !(
      metadata$sample_id %in%
         c(
            unlist(val_samples, use.names=F),
            rownames(pred_reports.cv$prob)
         )
   )
]

## Pred reports in one list
pred_reports <- c(
   list(CV=pred_reports.cv),
   lapply(val_samples, function(i){
      predict(
         object=model,
         newdata=features[i,],
         calc.feat.contrib=T, type='report', top.n.pred.classes=NULL
      )
   })
)

pred_reports <- lapply(pred_reports, function(i){
   i$class_actual <- as.factor(
      metadata$cancer_type[ match(rownames(i$prob), metadata$sample_id) ]
   )
   class(i) <- c('predReport',class(i))
   return(i)
})

# (function(){
#    male_samples <- rownames(features)[ features$gender.gender==FALSE ]
#    m <- pred_reports$CV$prob
#    m <- m[rownames(m) %in% male_samples,c('Ovary','Cervix','Uterus')]
#    m[order(-rowSums(m)),]
# })()
#
# (function(){
#    female_samples <- rownames(features)[ features$gender.gender==TRUE ]
#    m <- pred_reports$CV$prob
#    m <- m[rownames(m) %in% female_samples,'Prostate',drop=F]
#    m[order(-rowSums(m)),]
# })()


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
   df <- summary(pred_reports[[i]], top.n.classes=3, top.n.feat=5, prob.type='prob_scaled')
   data.frame(group=i, df)
}))

if(WRITE_OUTPUT){
   write.table(pred_summ, paste0(out_dir,'/pred_summ.txt'), sep='\t', quote=F, row.names=F)
}

## Low acc classes
# df <- calcFracCorrect(
#    actual=pred_reports$CV$class_actual,
#    predicted=pred_reports$CV$class_pred
# )
# dim(subset(df, frac_correct>=0.8))

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

## Patient report examples ================================
# View(subset(pred_summ, group=='holdout'))
if(WRITE_OUTPUT){
   ##
   anonymizeSampleId <- function(sample.id){
      #sample.id=rownames(features)
      cohort <- metadata$cohort[match(sample.id, metadata$sample_id)]
      is_hmf_sample <- cohort=='HMF'

      out <- sample.id
      out[is_hmf_sample] <- metadata$sample_id_2[match(sample.id[is_hmf_sample], metadata$sample_id)]

      return(out)
   }

   ##
   patientReportWrapper <- function(report, sample.name){
      #report=pred_reports$holdout
      #sample.name='XXXXXXXX'

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
      patientReportWrapper(pred_reports$holdout,'XXXXXXXX') ## Cervix
   )
   patient_reports <- cowplot::plot_grid(plotlist=patient_reports, ncol=1, rel_heights=c(1,1.6))

   pdf(paste0(out_dir,'/patient_reports.pdf'), 10, 11)
   patient_reports
   dev.off()

   # pdf(paste0(out_dir,'/patient_reports.pdf'), 10, 7)
   # patientReportWrapper(pred_reports$holdout,'DO35949') ## CNS_PiloAstro
   # patientReportWrapper(pred_reports$holdout,'XXXXXXXX') ## Skin_Other
   # dev.off()
}


## Main plots ================================
## Feature importance --------------------------------
if(WRITE_OUTPUT){

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

   ## --------------------------------
   p_imp_max <- maxImpPerFeatureType(pred_reports$CV$imp) +
      ylab('Cancer type') +
      theme(
         plot.margin=unit(c(0.01,0.04,0.01,0.01),'npc'),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         axis.text.x.top=element_blank(),
         axis.ticks.x.top=element_blank()
      )

   p_imp_sel <- (function(){
      sel_cancer_types <- c(
         'Cervix','CNS_PiloAstro','HeadAndNeck_Other','HeadAndNeck_SG',
         'Prostate','Sarcoma_GIST','Sarcoma_Lipo','Skin_Other'
      )
      m_imp <- model$imp[sort(sel_cancer_types),]

      topFeatures(m_imp, infer.feature.type=T, top.n=15, facet.ncol=2, feature.type.colors=FEATURE_TYPE_COLORS) +
         guides(fill=guide_legend(
            direction='horizontal',label.position='bottom',title.position='top',title.hjust=0.5,nrow=1,
            label.theme=element_text(angle=90, hjust=1, vjust=0.5, size=9)
         )) +
         theme(
            legend.position='bottom',
            legend.spacing.x=unit(0,'pt')
         )
   })()

   p_imp_combined <- cowplot::plot_grid(
      p_imp_max, p_imp_sel,
      align='h', axis='b', nrow=1, rel_widths=c(1, 0.8)
   )

   pdf(paste0(out_dir,'/imp_max.pdf'), 11, 8.5)
   plot(p_imp_combined)
   dev.off()

   ## --------------------------------
   p_imp_sel_list <- (function(){
      sel_cancer_types <- c('Cervix','Sarcoma_Lipo','Skin_Melanoma')
      m_imp <- model$imp[sort(sel_cancer_types),]
      topFeatures(
         m_imp, infer.feature.type=T, top.n=10, as.list=T, drop.legend.levels=F,
         feature.type.colors=FEATURE_TYPE_COLORS)
   })()

   pdf(paste0(out_dir,'/imp_list.pdf'), 4, 3)
   for(i in p_imp_sel_list) { plot(i) }
   dev.off()
}

## Confusion matrix --------------------------------
if(WRITE_OUTPUT){

   pdf(paste0(out_dir,'/confusion_heatmap.pdf'), 10.5, 10)

   plot(with(
      pred_reports$CV,
      confusionHeatmap2(
         actual=class_actual, probs=prob_scaled,
         rel.heights=c(100,11,22), plot.title='Cross-validation'
      )
   ))

   plot(with(
      pred_reports$holdout,
      confusionHeatmap2(
         actual=class_actual, probs=prob_scaled,
         rel.heights=c(100,11,22), plot.title='Holdout'
      )
   ))

   dev.off()

   # pdf(paste0(out_dir,'/confusion_heatmap.pdf'), 10.5, 10)
   #
   # plot(with(
   #    pred_reports$CV,
   #    confusionHeatmap(
   #       actual=class_actual, probs=prob_scaled,
   #       metrics=c('precision','recall'),
   #       which.plots=c('counts','perf','confusion'), rel.heights=c(0.31, 0.07, 1),
   #       plot.title='Cross-validation'
   #    )
   # ))
   #
   # plot(with(
   #    pred_reports$holdout,
   #    confusionHeatmap(
   #       actual=class_actual, probs=prob_scaled,
   #       metrics=c('precision','recall'),
   #       which.plots=c('counts','perf','confusion'), rel.heights=c(0.31, 0.07, 1),
   #       plot.title='Holdout'
   #    )
   # ))
   #
   # dev.off()
}

## Difference between CV and holdout perf --------------------------------
if(WRITE_OUTPUT){
   p_cv_holdout_perf_diff <- (function(){
      ## Calculate precision/recall
      calcPerfWrapper <- function(report){
         #report=pred_reports$CV

         require(mltoolkit)

         confusion <- confusionMatrix(
            predicted=oneHotEncode(report$class_pred),
            actual=report$class_actual,
            simplify=T
         )

         perf <- calcPerf(confusion, metrics=c('recall','precision'))
         perf$n_total <- confusion$tp + confusion$fn

         return(perf)
      }

      perf_holdout <- calcPerfWrapper(pred_reports$holdout)
      perf_cv <- calcPerfWrapper(pred_reports$CV)

      ## Calculate % difference
      percDiff <- function(x,y){ abs(x-y) / pmax(x,y) }
      perf_diff <- data.frame(
         class=perf_holdout$class,
         Recall=percDiff(perf_holdout$recall, perf_cv$recall),
         Precision=percDiff(perf_holdout$precision, perf_cv$precision)
      )
      perf_diff <- reshape2::melt(perf_diff, id.vars='class')
      colnames(perf_diff) <- c('class','metric','perc_diff')

      ## Calculate cohort totals
      n_total <- table(c(
         as.character(pred_reports$holdout$class_actual),
         as.character(pred_reports$CV$class_actual)
      ))
      n_holdout <- table(pred_reports$holdout$class_actual)

      perf_diff$n_total <- as.integer(n_total[ as.character(perf_diff$class) ])
      perf_diff$n_holdout <- as.integer(n_holdout[ as.character(perf_diff$class) ])

      ## Make label
      perf_diff$label <- paste0(perf_diff$class,' (',perf_diff$n_holdout,')')

      ## Plot
      ggplot(perf_diff, aes(x=n_total, y=perc_diff)) +
         facet_wrap(.~metric) +
         #xlim(-100,NA) +
         #ylim(-0.1,NA) +
         geom_point(color='red') +
         ggrepel::geom_text_repel(aes(label=label), size=2.5, segment.size=0.2, nudge_y=0.01, nudge_x=40) +
         ylab('% difference between\nCV and holdout performance') +
         xlab('x-axis: Total # of samples per cancer type\nLabel values: # of samples per cancer type in holdout set') +

         theme_bw() +
         theme(
            #panel.grid.minor=element_blank()
         )
   })()

   pdf(paste0(out_dir,'/cv_holdout_perf_diff.pdf'), 10, 7)
   plot(p_cv_holdout_perf_diff)
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
         'All' -> i_All$cancer_type_raw -> i_All$cancer_type_subgroup -> i_All$cancer_type_group
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
         i$color_index <- as.integer(droplevels(i$cancer_type_subgroup))
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
if(WRITE_OUTPUT){
   pdf(paste0(out_dir,'/topn_acc.pdf'), 11, 8)
   plot( with(pred_reports$CV, topnAcc(actual=class_actual, probs=prob_scaled)) + ggtitle('CV') )
   plot( with(pred_reports$holdout, topnAcc(actual=class_actual, probs=prob_scaled)) + ggtitle('holdout') )
   dev.off()
}

if(WRITE_OUTPUT){
   p_topn_acc <- (function(){
      df_cv <- with(pred_reports$CV, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))
      df_holdout <- with(pred_reports$holdout, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))

      df <- rbind(
         cbind(pred_type='CV', df_cv),
         cbind(pred_type='Holdout', df_holdout)
      )

      df$label <- with(df,paste0(
         round(frac,2),' (',count,'/',n_samples,')'
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
            name='Top-n class'
         ) +

         geom_text(
            aes(label=label, y=0.03),
            position=position_dodge(width=0.85),
            size=2.7, angle=90, hjust=0
         ) +

         scale_y_continuous(name='Top-n recall', limits=c(0,1)) +

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


## Plot RMD profiles ================================
if(WRITE_OUTPUT){
   p_rmd_profiles <- plotRmdProfile(t(model$rmd_sig_profiles))
   pdf(paste0(out_dir,'/rmd_profiles.pdf'), 11, 3)
   for(i in p_rmd_profiles){ plot(i) }
   dev.off()
}


# ## Validation performance ================================
# p_perf <- (function(){
#    df <- rbind(
#       with(pred_reports$CV,{
#          out <- calcFracCorrect(actual=class_actual, predicted=class_pred)
#          out$group <- 'Cross-validation'
#          return(out)
#       }),
#
#       with(pred_reports$holdout,{
#          out <- calcFracCorrect(actual=class_actual, predicted=class_pred)
#          out$group <- 'Indepedent test set'
#          return(out)
#       })
#    )
#
#    ## Class order
#    class_order <- with(subset(df,group=='Cross-validation'),{
#       class[ order(-frac_correct) ]
#    })
#    class_order <- unique(c('All',class_order))
#    df$class <- factor(df$class, class_order)
#    df <- df[order(df$class),]
#
#    ## Facet groups
#    facet.ncol=2
#    facet_n <- ceiling(length(class_order)/facet.ncol)
#    facet_groups <- rep(1:facet.ncol, each=facet_n)
#    facet_groups <- facet_groups[1:length(class_order)]
#    names(facet_groups) <- class_order
#    df$facet_group <- facet_groups[ as.character(df$class) ]
#
#    ## Label
#    df$label <- with(df,paste0(
#       round(frac_correct,2),' (',n_correct,'/',n_total,')'
#    ))
#    df$label <- factor(df$label, unique(df$label ))
#
#    ##
#    fill_colors <- structure(
#       c('#79AED2','#E1EBF6'),
#       names=levels(df$group)
#    )
#
#    ## Main
#    bar_width=0.8
#    ggplot(df, aes(x=class, y=frac_correct)) +
#       facet_wrap(.~facet_group, nrow=1, scales='free_y') +
#       geom_bar(
#          aes(fill=group), stat='identity',
#          position=position_dodge2(width=bar_width, preserve='single', reverse=T),
#          width=bar_width, color='black', size=0.25
#       ) +
#       geom_text(
#          aes(label=label, group=group), y=0.01, hjust=0, size=2.5,
#          position=position_dodge2(width=bar_width, preserve='single', reverse=T),
#       ) +
#
#       scale_fill_manual(values=fill_colors, name='Test group') +
#       #scale_fill_brewer(palette='Paired', name='Test group') +
#       scale_x_discrete(name='Cancer type', limits=rev) +
#       scale_y_continuous(name='Prop. correctly predicted', breaks=seq(0,1,0.2)) +
#       coord_flip() +
#       theme_bw() +
#       theme(
#          panel.grid.minor.x=element_blank(),
#          legend.position='top',
#          strip.background=element_blank(),
#          strip.text=element_blank()
#       )
# })()
#
# if(WRITE_OUTPUT){
#    pdf(paste0(out_dir,'/perf_cv_holdout.pdf'), 6, 6)
#    plot(p_perf)
#    dev.off()
# }
#
#
# ## Patient report ================================
# predReportSummaryCustom <- function(report){
#    #report <- pred_reports$holdout
#    summ <- predReportSummary(report)
#    summ$class_prob <- rowMaxs(report$prob_scaled)
#    summ$class_pred <- factor(
#       colnames(report$prob_scaled)[max.col(report$prob_scaled)],
#       levels=colnames(report$prob_scaled)
#    )
#
#    summ$is_correct_pred <- as.character(summ$class_actual) == as.character(summ$class_pred)
#    summ <- summ[order(summ$class_prob),]
#
#    return(summ)
# }
#
# predReportSummaryCustom(pred_reports$holdout)
# predReportSummaryCustom(pred_reports$CUP)
#
#
# p_patient_report <- (function(){
#    report <- pred_reports$CUP
#    #sample.name <- 'XXXXXXXX' ## Unknown, Uterus
#    #sample.name <- 'XXXXXXXX' ## Unknown, Gastric
#    sample.name <- 'XXXXXXXX' ## Unknown, Kidney
#    #sample.name <- 'XXXXXXXX'
#
#    plot.title <- paste0(
#       metadata$sample_id_2[ match(sample.name, metadata$sample_id) ], ## anon sample id
#       ' (actual cancer type: ',
#       metadata$cancer_type[ match(sample.name, metadata$sample_id) ],
#       ')'
#    )
#
#    p <- patientReport(
#       probs=report$prob_scaled[sample.name,],
#       feat.contrib=subset(report$feat_contrib, sample==sample.name),
#       plot.title=plot.title
#    )
#    plot(p)
#
#    pdf(paste0(model_dir,'/plots/patient_report_example.pdf'), 6, 5.5)
#    grid::grid.draw(p)
#    dev.off()
# })()
#
#
#
#
