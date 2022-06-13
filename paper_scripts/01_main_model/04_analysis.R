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
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))
#library(matrixStats)

## Load data ================================
## Model --------------------------------
model_dir <- paste0(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/'))
out_dir <- paste0(model_dir,'/report/')
dir.create(out_dir, showWarnings=F)

metadata <- read.delim(paste0(model_dir,'/features/metadata_training.txt'))

features <- cacheAndReadData(paste0(model_dir,'/features/features_allSamples.txt.gz'), overwrite=T)

# (function(){
#    devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/featureExtractor/'))
#    rmd_profiles <- read.delim(RMD_BINS_HG19_PATH)
#    m <- fitToRmdProfiles(features, rmd.profiles=model$rmd_sig_profiles)
#    m_ss <- m[
#       subset(metadata, is_training_sample | is_holdout_sample, sample_id, drop=T),
#       grep('rmd',colnames(m))
#     ]
#    write.table(m_ss, gzfile('~/Desktop/rmd_sig_contribs.txt.gz'), sep='\t', quote=F)
# })()

anonymizeSampleId <- function(sample.id){
   #sample.id=pred_summ$sample
   #cohort <- metadata$cohort[match(sample.id, metadata$sample_id)]
   #is_hmf_sample <- cohort=='HMF'
   is_hmf_sample <- grepl('^(CPCT|WIDE|DRUP|ACTN)',sample.id)

   out <- sample.id
   out[is_hmf_sample] <- metadata$sample_id_2[match(sample.id[is_hmf_sample], metadata$sample_id)]

   return(out)
}

model <- cacheAndReadData(paste0(model_dir,'/final/model.rds'), overwrite=T)

## Pred reports --------------------------------
## CV set
gatherCvOutput(
   cv.out.dir=paste0(model_dir,'/cv_out/'),
   out.path=paste0(model_dir,'/test_set_report.rds')
)
pred_reports.cv <- cacheAndReadData(paste0(model_dir,'/test_set_report.rds'), overwrite=T)

## Val set samplesLuan Nguyen <N.L.Nguyen-2@umcutrecht.nl>
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
      out <- predict.randomForestEnsemble(
         object=model,
         newdata=features[i,],
         calc.feat.contrib=T,
         type='report', top.n.pred.classes=NULL
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

## Main perf plots #################################################################################
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

if(WRITE_OUTPUT){
   write.table(pred_summ, paste0(out_dir,'/pred_summ.txt'), sep='\t', quote=F, row.names=F)
}

## Prob calibration --------------------------------
## Calibration curves
if(WRITE_OUTPUT){
   p_cal_curves <- probCalCurves(report=pred_reports$CV, method='isotonic', prob.prescale=T, output='plot', facet.ncol=5)
   pdf(paste0(out_dir,'/prob_calib_curves.pdf'), 11, 11)
   plot(p_cal_curves)
   dev.off()

   write.table(cal_curves, paste0(out_dir,'/prob_calib_curves.txt'), sep='\t', quote=F, row.names=F)
}

## Reliability curves
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
      df2$class <- paste0(df2$class,' ')

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
## Performance --------------------------------
if(WRITE_OUTPUT){

   ## Init --------------------------------
   df <- subset(
      pred_summ,
      group %in% c('CV','holdout'),
      c(group, sample, actual_class, pred_correct, pred_class.1, pred_class.2, pred_class.3)
   )
   df$pred_class.1 <- as.factor(df$pred_class.1)
   df$actual_class <- droplevels(df$actual_class)

   ## Sample size --------------------------------
   ## Count per cancer type
   pd.samples <- as.data.frame(table(actual=df$actual_class, group=df$group))

   ## Count for all samples per cohort
   pd.samples <- rbind(
      data.frame(
         actual='All',
         group=c('CV','holdout'),
         Freq=as.integer(table(df$group))
      ),
      pd.samples
   )

   ## Plot
   ymax <- 500
   bar_width <- 0.8

   vlines <- unique(as.integer(as.factor(pd.samples$actual))) - 0.5
   vlines <- vlines[-1]
   hlines <- unique(as.integer(as.factor(pd.samples$actual))) - 0.5

   p.samples <- ggplot(pd.samples, aes(x=actual, y=Freq, fill=group)) +

      geom_vline(xintercept=vlines, color='grey50', size=0.3) +
      geom_bar(
         stat='identity', position=position_dodge(width=bar_width), width=bar_width,
         color='grey50', size=0.2
      ) +
      geom_text(
         aes(label=Freq, y=0.02*ymax), position=position_dodge(width=bar_width),
         angle=90, hjust=0, vjust=0.5, size=2.5
      ) +

      scale_fill_manual(
         values=c(CV='#9FCADF',holdout='#4693C3'),
         labels=c(CV='Train set',holdout='Test set'),
         guide=guide_legend(override.aes=list(color=NA))
      ) +

      scale_x_discrete(position='top', name='Actual class', expand=c(0.014, 0.014)) +
      scale_y_continuous(name='# samples') +
      coord_cartesian(ylim=c(0, ymax)) +

      theme_bw() +
      theme(
         plot.margin=unit(c(0, 5, 5, 5),'pt'),
         panel.grid=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         legend.position='left',
         #legend.direction='horizontal',
         legend.title=element_blank()
      )



   ## Performance summary --------------------------------
   ## Functions
   calcPrecision <- function(df, group.name=NULL){
      ## Calculate no. correct for each predicted class
      l <- lapply(levels(df$pred_class.1), function(i){
         #i='Breast'
         i_df <- df[df$pred_class.1==i,]
         n_correct <- sum(i_df$actual_class==i_df$pred_class.1)
         n_total <- nrow(i_df)

         data.frame(
            class=i,
            n_correct,
            n_total
         )
      })
      perf <- do.call(rbind, l)

      ## Add pan-cancer performance
      perf <- rbind(
         data.frame(
            class='All',
            n_correct=sum(perf$n_correct),
            n_total=sum(perf$n_total)
         ),
         perf
      )

      ## Calculate precision
      perf$frac_correct <- perf$n_correct/perf$n_total
      perf$metric <- 'Precision'

      return(perf)
   }

   calcTopNRecall <- function(df, top.n=1){
      agg <- aggregate(
         df$pred_correct<=top.n,
         list(class=df$actual_class),
         function(x){ c(n_correct=sum(x), n_total=length(x)) }
      )

      perf <- data.frame(class=agg$class, agg$x)

      ## Add pan-cancer performance
      perf <- rbind(
         data.frame(
            class='All',
            n_correct=sum(perf$n_correct),
            n_total=sum(perf$n_total)
         ),
         perf
      )

      ## Calculate recall
      perf$frac_correct <- perf$n_correct/perf$n_total
      if(top.n==1){
         perf$metric <- 'Recall'
      } else {
         perf$metric <- sprintf('Top-%s recall', top.n)
      }

      return(perf)
   }

   ## Prep plot data
   pd.perf_summ <- do.call(rbind, lapply(unique(df$group), function(i){
      i_df <- subset(df, group==i)
      out <- rbind(
         calcTopNRecall(i_df, top.n=1),
         calcTopNRecall(i_df, top.n=2),
         calcPrecision(i_df)
      )
      cbind(group=i, out)
   }))

   pd.perf_summ$perc_correct <- round(pd.perf_summ$frac_correct*100)
   #pd.perf_summ$label <- paste0(pd.perf_summ$perc_correct,' (',pd.perf_summ$n_correct,')')
   #pd.perf_summ$label <- paste0(pd.perf_summ$perc_correct,'%')
   pd.perf_summ$metric <- factor(pd.perf_summ$metric, unique(pd.perf_summ$metric))

   p.perf_summ <- ggplot(pd.perf_summ, aes(x=group, y=metric, fill=frac_correct, group=group)) +
      facet_grid(~class) +
      geom_tile() +
      geom_text(aes(label=perc_correct), angle=90, size=2.5) +
      geom_hline(yintercept=hlines, color='black', size=0.3) +
      scale_fill_distiller(palette='RdYlGn', direction=1) +
      scale_x_discrete(name='Group', expand=c(0,0)) +
      scale_y_discrete(name='Metric', expand=c(0,0), limits=rev) +
      theme_bw() +
      theme(
         plot.margin=unit(c(0, 5, 5, 5),'pt'),
         panel.grid=element_blank(),
         panel.spacing=unit(0, 'pt'),
         panel.border=element_rect(color='black', size=0.3),

         #strip.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         strip.text.x=element_blank(),
         strip.background.x=element_blank(),
         axis.title=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank(),

         legend.position='none'
      )


   ## Confusion matrix --------------------------------
   confusionMatrixPlotData <- function(actual, predicted, rel.values=T, output.format='long'){
      # df_ss <- subset(pred_summ, group=='CV')
      # actual=df_ss$actual_class
      # predicted=df_ss$pred_class.1

      ## Init ----------------------------
      if(!(output.format %in% c('wide','long'))){
         stop("`output.format` must be 'wide' or 'long'")
      }

      predicted_classes <- sort(unique(as.character(predicted)))

      if(!is.factor(actual)){ stop('`actual` must be a factor') }

      ## Only show classes that are present in predicted
      is_predicted_class <- actual %in% predicted_classes
      actual <- actual[is_predicted_class]
      predicted <- predicted[is_predicted_class]
      rm(is_predicted_class)

      ##
      actual <- factor(actual, predicted_classes)
      predicted <- factor(predicted, predicted_classes)

      ## Large confusion matrix ----------------------------
      cm_raw <- table(predicted, actual)
      actual_counts <- table(actual)

      ## Only keep predicted classes on prediction axis
      cm <- cm_raw
      cm <- cm[predicted_classes,]

      ## Convert to fraction classified as which class
      overall_acc <- sum( as.character(predicted)==as.character(actual) )
      if(rel.values){
         cm <- apply(cm,2,function(i){ i/sum(i) })
         cm <- round(cm,2)
         overall_acc <- round(overall_acc / sum(actual_counts), 2)
      }

      ## Add overall accuracy
      cm <- cbind(All=NaN,cm)
      cm <- rbind(All=NaN,cm)
      cm[1,1] <- overall_acc

      if(output.format=='wide'){ return(cm) }

      out <- reshape2::melt(cm)
      colnames(out) <- c('predicted','actual','value')

      return(out)
   }

   ## Plot data
   pd.cm <- lapply(c('CV','holdout'), function(i){
      df_ss <- subset(df, group==i)
      out <- confusionMatrixPlotData(df_ss$actual_class, df_ss$pred_class.1)
      out <- cbind(group=i, out)
      return(out)
   })
   pd.cm <- do.call(rbind, pd.cm)
   pd.cm$group <- factor(pd.cm$group, unique(pd.cm$group))

   ## Plot
   p.cm <- ggplot(pd.cm, aes(x=group, y=predicted)) +
      facet_grid(~actual, switch='x') +

      geom_tile(aes(fill=value), width=1.05) +
      geom_text(aes(label=round(value*100)), angle=90, size=2.5) +
      geom_hline(yintercept=hlines, color='grey50', size=0.3) +

      scale_fill_distiller(palette='YlGnBu', na.value='grey75', direction=-1) +
      scale_x_discrete(name='Actual class', expand=c(0,0), position='bottom') +
      scale_y_discrete(name='Predicted class', expand=c(0,0), limits=rev) +

      theme_bw() +
      theme(
         plot.margin=unit(c(0, 5, 5, 5),'pt'),
         panel.grid=element_blank(),
         panel.spacing=unit(0, 'pt'),
         panel.border=element_rect(color='grey50', size=0.3),

         strip.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         strip.background.x=element_blank(),

         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank(),

         legend.position='none'
      )


   ## Combine ----------------------------
   # p.combined <- cowplot::plot_grid(
   #    p.samples, p.perf_summ, p.cm,
   #    ncol=1, rel_heights=c(30, 23, 100), align='v', axis='tblr'
   # )

   p.combined <- cowplot::plot_grid(
      p.samples, p.perf_summ, p.cm,
      ncol=1, rel_heights=c(30, 8, 100), align='v', axis='tblr'
   )

   pdf(paste0(out_dir,'/confusion_heatmap.split.pdf'), 12, 13)
   plot(p.combined)
   dev.off()

   ## Write raw data ----------------------------
   confusion_heatmap_values <- lapply(c('CV','holdout'), function(i){
      #i='CV'
      df_ss <- subset(df, group==i)
      out <- confusionMatrixPlotData(df_ss$actual_class, df_ss$pred_class.1, output.format='wide', rel.values=F)
      out <- data.frame('predicted class \\ actual class'=rownames(out), out, check.names=F)
      out[is.na(out)] <- NA
      rownames(out) <- NULL
      return(out)
   })
   names(confusion_heatmap_values) <- c('confusion_matrix.CV', 'confusion_matrix.holdout')
   confusion_heatmap_values$performance_metrics <- subset(pd.perf_summ, select=-perc_correct)
   openxlsx::write.xlsx(confusion_heatmap_values, paste0(out_dir,'/confusion_heatmap.values.xlsx'))
}

if(WRITE_OUTPUT){
   perfCurveData <- function(report=NULL, prob=NULL, actual=NULL, metrics=c('ppv','tpr')){
      if(F){
         prob=pred_reports$CV$prob_scaled
         actual=pred_reports$CV$class_actual
         metrics=c('tpr','ppv')
      }

      ##
      if(!is.null(report)){
         prob <- report$prob_scaled
         actual <- report$class_actual
      }

      ## Calc perf
      confusion <- mltoolkit::confusionMatrix(predicted=prob, actual=actual, cutoff.interval=0.01)
      perf <- do.call(rbind, lapply(names(confusion), function(i){
         df <- mltoolkit::calcPerf(confusion=confusion[[i]], metrics=metrics, add.start.end.values=F)
         df$class <- i
         df$n_samples <- nrow(df)
         return(df)
      }))

      ## Plot data
      pd <- reshape2::melt(perf, measure.vars=metrics)
      colnames(pd)[colnames(pd)=='variable'] <- 'metric'
      pd$metric <- toupper(pd$metric)
      pd$metric <- factor(pd$metric, unique(pd$metric))
      return(pd)
   }

   pd <- (function(){
      pd1 <- perfCurveData(pred_reports$CV)
      pd1$group <- 'CV'
      pd2 <- perfCurveData(pred_reports$holdout)
      pd2$group <- 'holdout'

      rbind(pd1, pd2)
   })()

   pd <- pd[order(pd$metric, pd$group),]
   pd$metric_group <- paste0(pd$metric,' | ', pd$group)
   pd$metric_group <- factor(pd$metric_group, unique(pd$metric_group))

   p <- ggplot(pd, aes(x=cutoff, y=value, color=metric_group, linetype=metric_group)) +
      facet_wrap(~class) +
      geom_path(size=0.4) +
      scale_color_manual(name='Metric | Group', values=c('#984EA3','#984EA3','#FF7F00','#FF7F00')) +
      scale_linetype_manual(name='Metric | Group', values=c('solid','11','solid','11')) +
      guides(color=guide_legend(override.aes=list(size=0.7))) +
      labs(y='Metric value', x='Probability cutoff') +
      theme_bw() +
      theme(
         panel.grid.minor=element_blank()
      )

   pdf(paste0(out_dir,'/perf_curves.pdf'), 14, 10)
   plot(p)
   dev.off()

}

## Feature importance --------------------------------
if(WRITE_OUTPUT){

   ## Plots
   plots <- list()
   plots[[1]] <- topFeatures(model$imp, infer.feature.type=T, top.n=15, facet.ncol=5, feature.type.colors=FEATURE_TYPE_COLORS) + ggtitle('Final model')
   plots[[2]] <- topFeatures(pred_reports$CV$imp, infer.feature.type=T, top.n=15, facet.ncol=5, feature.type.colors=FEATURE_TYPE_COLORS) + ggtitle('CV')

   pdf(paste0(out_dir,'/imp_top.pdf'), 12, 12)
   for(i in plots){ plot(i) }
   dev.off()

   ## Export table
   write.table(t(model$imp), paste0(out_dir,'imp.txt'), sep='\t', quote=F)
}

## Misc ############################################################################################
## Sample counts ================================
if(WRITE_OUTPUT){
   metadata_ss <- subset(metadata, is_holdout_sample | is_training_sample)
   metadata_ss$cancer_type <- as.factor(metadata_ss$cancer_type)

   ## Count train/test samples
   cohort_names <- c('All','HMF','PCAWG')
   counts <- do.call(rbind, lapply(cohort_names, function(i){
      #i='All'
      if(i=='All'){
         df <- metadata_ss
      } else {
         df <- subset(metadata_ss, cohort==i)
      }

      tab <- with(
         df,
         table(cancer_type=cancer_type, is_training_sample=is_training_sample)
      )
      tab <- as.data.frame(unclass(tab))
      colnames(tab) <- c('n_test','n_train')
      tab <- data.frame(
         cohort=i,
         cancer_type=rownames(tab),
         n_total=tab$n_test + tab$n_train,
         n_train=tab$n_train,
         n_test=tab$n_test
      )
      rownames(tab) <- NULL

      return(tab)
   }))
   counts$frac_test <- counts$n_test/counts$n_total
   #counts$frac_test[is.na(counts$frac_test)] <- 0

   ## Format table for exporting
   fixed_cols <- c('cohort','cancer_type')
   counts_export <- lapply(cohort_names, function(i){
      #i='All'
      df <- subset(counts, cohort==i)
      #df_fixed <- df[,fixed_cols]
      df_values <- df[,!(colnames(df) %in% fixed_cols)]
      df_values$frac_test <- round(df_values$frac_test, 3)
      colnames(df_values) <- paste0(colnames(df_values),'.',i)
      return(df_values)
   })
   counts_export <- do.call(cbind, counts_export)

   counts_export <- cbind(
      counts[counts$cohort=='All',fixed_cols],
      counts_export
   )
   counts_export$cohort <- NULL

   ## Write output
   write.table(
      counts_export, paste0(out_dir,'/sample_counts.txt'),
      sep='\t',quote=F,row.names=F
   )
}

## Patient report examples ================================
# View(subset(pred_summ, group=='holdout'))
if(WRITE_OUTPUT){

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
            drop.feature.type.levels=F
         )
      })
   }

   # with(pred_reports$holdout,{
   #    patientReport(
   #       probs=prob_scaled,
   #       feat.contrib=feat_contrib,
   #       sample.name='XXXXXXXX',
   #       #sample.name='DO36039',
   #       #plot.title=plot_title,
   #       prob.thres.min=1,
   #       rel.widths=c(1, 1, 1), drop.feature.type.levels=F
   #    )
   # })
   #View(subset(pred_summ, group %in% c('CV','holdout')))
   #View( subset(pred_summ, group %in% c('holdout') & grepl('HMF',sample_anon) & pred_class.1>=0.7) )
   patientReportWrapper(pred_reports$holdout, sample.name='XXXXXXXX')

   patient_reports <- list(
      #patientReportWrapper(pred_reports$holdout, sample.name='DO35997'), ## CNS_PiloAstro
      patientReportWrapper(pred_reports$holdout, sample.name='XXXXXXXX'), ## Lung_NonSmallCell
      patientReportWrapper(pred_reports$CV,'DO7304') ## Lymphoid

   )
   patient_reports <- cowplot::plot_grid(plotlist=patient_reports, ncol=1, rel_heights=c(1,1.5))

   pdf(paste0(out_dir,'/patient_reports.pdf'), 11, 11)
   patient_reports
   dev.off()


   if(F){
      ## For advanced omics course
      sel_sample_ids <- c('DO220857','DO25189','XXXXXXXX','DO28089','DO42400','XXXXXXXX')
      lp <- lapply(sel_sample_ids, function(i){
         patientReportWrapper(pred_reports$CV, i)
      })
      pdf(paste0(out_dir,'/patient_reports.adv_omics_course.pdf'), 11, 5)
      for(i in lp){ plot(i) }
      dev.off()
   }
}

## Features ================================
## Feature importance of selected classes --------------------------------
if(WRITE_OUTPUT){

   # imp_exp <- t(model$imp)
   # imp_exp <- data.frame(feature=rownames(imp_exp), imp_exp, row.names=NULL)

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
         'Cervix','CNS_PiloAstro','HeadAndNeck_Other','Liver',
         'Lung_NonSmallCell','Prostate','Sarcoma_Lipo','Skin_Carcinoma'
      )
      m_imp <- model$imp[sort(sel_cancer_types),]

      topFeatures(m_imp, infer.feature.type=T, top.n=15, facet.ncol=2, feature.type.colors=FEATURE_TYPE_COLORS) +
         guides(fill=guide_legend(
            direction='horizontal',label.position='bottom',title.position='top',title.hjust=0.5,nrow=1,
            label.theme=element_text(angle=90, hjust=1, vjust=0.5, size=9)
         )) +
         theme(
            #plot.margin=unit(c(0.01,0.04,0.01,0.01),'npc'),
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

## Feature averages per class --------------------------------
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
      feature_names <- insValueBefore(feature_names,'rmd.Pancreas.1','rmd.Pancreas.2')
      feature_names <- feature_names[ -head(which(feature_names=='rmd.Pancreas.1'),1) ]

      feature_names <- insValueBefore(feature_names,'rmd.Sarcoma_Other.1','rmd.Skin_Carcinoma.1')
      feature_names <- feature_names[ -head(which(feature_names=='rmd.Sarcoma_Other.1'),1) ]

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
               panel.spacing = unit(2,'pt'),
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

## Performance ================================
## Raw vs scaled probs --------------------------------
if(WRITE_OUTPUT){
   p_raw_vs_scaled_probs <- (function(){
      l <- list(
         CV.raw = pred_reports$CV[c('class_actual','prob')],
         CV.calibrated = pred_reports$CV[c('class_actual','prob_scaled')],
         holdout.raw = pred_reports$holdout[c('class_actual','prob')],
         holdout.calibrated = pred_reports$holdout[c('class_actual','prob_scaled')]
      )

      df <- do.call(rbind, lapply(names(l), function(i){
         #i='CV.raw'
         probs <- l[[i]][[2]]
         class_predicted <- colnames(probs)[max.col(probs)]
         class_predicted <- factor(class_predicted, colnames(probs))
         out <- calcFracCorrect(actual=l[[i]]$class_actual, predicted=class_predicted)
         cbind(group=i, out)
      }))
      df$val_type <- sub('[.].+$','',df$group)
      df$prob_type <- sub('^.+[.]','',df$group)
      df$prob_type <- factor(df$prob_type, unique(df$prob_type))

      df$label <- paste0(
         round(df$frac_correct, 2),
         ' (',df$n_correct,'/',df$n_total,')'
      )

      p <- ggplot(df, aes(x=val_type, y=frac_correct, fill=prob_type)) +
         facet_wrap(class~.) +
         geom_bar(stat='identity', position=position_dodge(width=0.9), color='black', size=0.3) +
         geom_text(
            aes(label=label, y=0.02), position=position_dodge(width=0.9),
            angle=90, hjust=0, vjust=0.5, size=2.7
         ) +
         labs(y='Recall', x='Validation type', fill='Prob. type') +
         scale_fill_manual(values=c(raw='#b7d7e8', calibrated='#4b93c3')) +
         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            legend.position='bottom'
         )
   })()

   pdf(paste0(out_dir,'/raw_vs_cal_probs.pdf'), 11, 10)
   plot(p_raw_vs_scaled_probs)
   dev.off()
}

## Difference between CV and holdout perf --------------------------------
if(WRITE_OUTPUT){

   calcRecallAndPrecision <- function(pred_report){
      #pred_report <- pred_reports$CV

      df <- data.frame(
         actual=pred_report$class_actual,
         predicted=pred_report$class_pred
      )

      ## Recall
      recall <- aggregate(
         df$actual==df$predicted,
         list(class=df$actual),
         function(x) c(n_correct=sum(x), n_total=length(x))
      )
      recall <- as.data.frame(recall)
      recall <- data.frame(metric='Recall', class=recall$class, recall$x)

      ## Precision
      precision <- aggregate(
         df$actual==df$predicted,
         list(class=df$predicted),
         function(x) c(n_correct=sum(x), n_total=length(x))
      )
      precision <- as.data.frame(precision)
      precision <- data.frame(metric='Precision', class=precision$class, precision$x)

      ##
      out <- rbind(recall, precision)
      out$frac_correct <- out$n_correct / out$n_total
      return(out)
   }

   p_diff_cv_holdout <- (function(){
      ## Calc perf
      df_cv <- calcRecallAndPrecision(pred_reports$CV)
      df_cv$group <- 'CV'

      df_h <- calcRecallAndPrecision(pred_reports$holdout)
      df_h$group <- 'holdout'

      ## Merge CV and holdout perf
      df <- data.frame(
         class=df_cv$class,
         metric=df_cv$metric,

         frac_correct.cv=df_cv$frac_correct,
         frac_correct.holdout=df_h$frac_correct,
         n_total.cv=df_cv$n_total,
         n_total.holdout=df_h$n_total
      )
      df$metric <- factor(df$metric, unique(df$metric))
      df$n_total <- df$n_total.cv + df$n_total.holdout

      ## Calc diff between holdout and CV accuracy
      relDiff <- function(x,y){ abs(x-y) / pmax(x,y) }
      df$frac_correct.perc_diff <- relDiff(df$frac_correct.h, df$frac_correct.cv)

      ## Point labels
      df$label <- paste0(df$class,' (',df$n_total.holdout,')')

      ## Plot
      ggplot(df, aes(x=n_total, y=frac_correct.perc_diff)) +
         facet_grid(~metric) +
         geom_point(color='red') +
         ggrepel::geom_text_repel(aes(label=label), size=2.5, segment.size=0.2, nudge_y=0.01, nudge_x=40) +
         ylab('Relative difference in performance\nbetween CV and holdout sets') +
         xlab('x-axis: Total no. of samples per cancer type\nLabel values: No. of samples per cancer type in holdout set') +
         theme_bw()  +
         theme(
            panel.grid.minor=element_blank()
         )
   })()

   pdf(paste0(out_dir,'/diff_cv_holdout.pdf'), 10, 7)
   plot(p_diff_cv_holdout)
   dev.off()




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
         ylab('Relative difference in recall\nbetween CV and holdout sets') +
         xlab('x-axis: Total no. of samples per cancer type\nLabel values: No. of samples per cancer type in holdout set') +
         theme_bw()
   })()

   pdf(paste0(out_dir,'/diff_cv_holdout.pdf'), 7, 6)
   plot(p_diff_cv_holdout)
   dev.off()
}

# ## Top-N accuracy, CV and holdout sets --------------------------------
# if(WRITE_OUTPUT){
#    p_topn_acc <- (function(){
#       df_cv <- with(pred_reports$CV, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))
#       df_holdout <- with(pred_reports$holdout, topnAcc(actual=class_actual, probs=prob_scaled, output='values'))
#
#       df <- rbind(
#          cbind(pred_type='CV', df_cv),
#          cbind(pred_type='Holdout', df_holdout)
#       )
#
#       df$label <- with(df,paste0(
#          round(frac,2),' (',correct,'/',total,')'
#       ))
#
#       ggplot(df, aes(x=pred_type, y=frac, group=class_num)) +
#          facet_wrap(~actual) +
#          geom_bar(
#             aes(fill=class_num),
#             stat='identity', position='dodge',
#             color='black', size=0.3, width=0.85
#          ) +
#          scale_fill_gradient(
#             low='#4390BC', high='#EFF6B9', guide='legend', breaks=unique(df$class_num),
#             name='Top-N class'
#          ) +
#
#          geom_text(
#             aes(label=label, y=0.03),
#             position=position_dodge(width=0.85),
#             size=2.7, angle=90, hjust=0
#          ) +
#
#          scale_y_continuous(name='Top-N recall', limits=c(0,1)) +
#          xlab('Validation type') +
#
#          theme_bw() +
#          theme(
#             panel.grid.minor=element_blank(),
#             panel.grid.major.x=element_blank(),
#             legend.position='bottom'
#          )
#    })()
#
#    pdf(paste0(out_dir,'/topn_acc.pdf'), 11, 10)
#    plot(p_topn_acc)
#    dev.off()
# }

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
         #perf$n_samples <- confusion$tp + confusion$fn

         cbind(
            confusion,
            n_samples=confusion$tp + confusion$fn,
            perf[,metric_names,drop=F]
         )
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
      i_perf <- subset(perf_comparison, metric==i & !is.na(metric_value))
      ggplot(i_perf, aes(x=classifier_name, y=metric_value, group=cancer_type_subgroup)) +
         #facet_grid(classifier_name~cancer_type_group, scales='free_x') +
         facet_wrap(cancer_type_group~., ncol=5) +

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

   #lp_perf_comparison$recall <- lp_perf_comparison$recall + ylab('Fraction of samples correctly predicted')

   pdf(paste0(out_dir,'/other_classifiers_comparison.pdf'), 12, 12)
   for(i in lp_perf_comparison){ plot(i) }
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

      ## Add MSI, HRD status, etc
      metadata_full <- read.delim(paste0(model_dir,'/features/metadata_full.txt'))
      df <- cbind(
         df,
         metadata_full[
            match(df$sample, metadata_full$sample_id),
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
}

## Impact of coverage ================================
if(WRITE_OUTPUT){
   ##
   metadata_full <- read.delim(paste0(model_dir,'/features/metadata_full.txt'), na.strings=c('NA','#N/A'))

   pred_summ_ss <- subset(
      pred_summ,
      group %in% c('CV','holdout'),
      #group=='CV',
      c(group,sample,sample_anon,actual_class,pred_correct,pred_class.1,pred_class.2,pred_class.3,pred_prob.1,pred_prob.2,pred_prob.3)
   )
   pred_summ_ss$coverage_bin <- metadata_full$coverage_bin[ match(pred_summ_ss$sample, metadata_full$sample_id) ]

   ##
   calcFracCorrectByCoverage <- function(df){

      df <- rbind(
         within(df, actual_class <- 'All'),
         df
      )

      agg <- with(df,{
         aggregate(
            pred_correct==1,
            list(actual_class=actual_class, coverage_bin=coverage_bin),
            function(x){ c(n_total=length(x), n_correct=sum(x), frac_correct=sum(x)/length(x)) }
         )
      })
      agg <- cbind(agg[-length(agg)], agg$x)
      agg$label <- paste0(round(agg$frac_correct,2),' (',agg$n_correct,'/',agg$n_total,')')

      agg$sample_size_bin <- cut(
         agg$n_total,
         breaks=c(0,10,20,Inf), labels=c('0-10','11-20','>=21'),
         include.lowest=T
      )

      return(agg)
   }

   pd_bars <- calcFracCorrectByCoverage(subset(pred_summ_ss, group=='CV'))

   sel_ct <- lapply(split(pd_bars$n_total, pd_bars$actual_class), function(i){ sum(i>=3)==3 })
   sel_ct <- names(sel_ct)[sel_ct==TRUE]

   p <- ggplot(subset(pd_bars, actual_class %in% sel_ct), aes(x=coverage_bin, y=frac_correct)) +
      facet_wrap(~actual_class) +

      geom_bar(aes(fill=coverage_bin), stat='identity', size=0.3, color='black') +
      scale_fill_brewer(palette='Blues') +

      geom_text(
         aes(label=label, y=0.02),
         angle=90, hjust=0, vjust=0.5, size=2.5
      ) +

      labs(
         y='Recall',
         fill='Tumor WGS\ncoverage bin'
      ) +

      theme_bw() +
      theme(
         panel.grid.minor=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank()
      )

   pdf(paste0(out_dir,'/coverage_comparison.pdf'), 9, 7)
   plot(p)
   dev.off()
}


## Write features ================================
if(WRITE_OUTPUT){
   m <- features
   rownames(m) <- anonymizeSampleId(rownames(m))
   write.table(m, gzfile('/Users/lnguyen/Desktop/HMF_PCAWG.cuplr_features.txt.gz'), sep='\t', quote=F)
}

## Plot RMD profiles ================================
if(WRITE_OUTPUT){
   m <- model$rmd_sig_profiles
   m <- m[,order(colnames(m))]

   p_rmd_profiles <- plotRmdProfile(t(m))
   pdf(paste0(out_dir,'/rmd_profiles.pdf'), 11, 3)
   for(i in p_rmd_profiles){ plot(i) }
   dev.off()
}

## Confidence score ######################################################################################
multiclassCombinedScore <- function(probs, top.n=NULL){
   if(F){
      #probs=pred_reports$CV$prob_scaled
      probs=pred_reports$holdout$prob_scaled['XXXXXXXX',]
      top.n=3
   }


   ## Checks --------------------------------
   if(is.vector(probs)){
      probs <- matrix(probs, nrow=1, dimnames=list(NULL, names(probs)))
   }

   ## Main --------------------------------
   ## Per sample; sort probs
   ##
   ## 1.0 |
   ##     |o - - - - - - - - -
   ## 0.8 |o|                 | <- box = non_prob1_total_area
   ##     |o|                 | <- empty space in box = non_prob1_non_bar_area
   ## 0.6 |o|                 |    -> Then normalize between 0-1: final certainty score
   ##     |o|                 |
   ## 0.4 |o|                 |
   ##     |o|o                |
   ## 0.2 |o|o                |
   ##     |o|o o              |
   ## 0   |o|o o o            |
   ##      --------------------
   ##      1 2 3 4 5 ...
   ##
   ##

   #t(t(sort(probs['XXXXXXXX',], decreasing=T)))

   probs_descending <- t(apply(probs,1,sort,decreasing=T))

   if(!is.null(top.n)){
      probs_descending <- probs_descending[,1:top.n,drop=F]
   }

   prob_1 <- probs_descending[,1]
   n_classes_minus_1 <- ncol(probs_descending) - 1

   rowSums(prob_1 - probs_descending[,-1,drop=F]) / n_classes_minus_1

   # ####
   # probs_descending <- t(apply(probs,1,sort,decreasing=T))
   #
   # prob_1 <- probs_descending[,1]
   # n_classes_minus_1 <- ncol(probs) - 1
   #
   # non_prob1_total_area <- prob_1 * n_classes_minus_1
   # non_prob1_bar_area <- rowSums(probs_descending[,-1,drop=F])
   # non_prob1_non_bar_area <- non_prob1_total_area - non_prob1_bar_area
   #
   # certainty <- non_prob1_non_bar_area / n_classes_minus_1
   # return(certainty)
}

function(report=NULL, actual=NULL, probs=NULL){
   if(F){
      report <- pred_reports$CV

      actual=report$class_actual
      probs=report$prob_scaled
   }

   # sum(
   #    report$class_actual==colnames(report$prob_scaled)[max.col(report$prob_scaled)]
   # )

   ## Init --------------------------------
   if(!is.null(report)){
      actual <- report$class_actual

      if(!is.null(report$prob_scaled)){
         probs <- report$prob_scaled
      } else {
         probs <- report$prob
      }

   }

   if(length(actual)!=nrow(probs)){
      stop('length(actual) and nrow(probs) must be equal')
   }

   if(!is.factor(actual)){
      warning('`actual` is not a factor. Factor levels will be set automatically')
      actual <- as.factor(actual)
   }
   names(actual) <- rownames(probs)

   ## Main --------------------------------
   uniq_classes <- colnames(probs)
   top_classes <- t(apply(probs,1,function(i){
      uniq_classes[order(i, decreasing=T)]
   }))

   top_correct <- sweep(top_classes, 1, as.character(actual), '==')

   topn_correct <- t(apply(top_correct, 1, function(i){
      which_true <- which(i)
      if(length(which_true)!=0){
         i[ which_true:length(i) ] <- TRUE
      }
      return(i)
   }))

   ## --------------------------------
   mccs <- multiclassCombinedScore(probs)
   #mccs <- probs[,1]


   # multiclassCombinedScore(
   #    pred_reports$holdout$prob_scaled['XXXXXXXX',]
   # )
   #
   # library(magrittr)
   #
   # pred_reports$holdout$prob_scaled['XXXXXXXX',] %>%
   #    sort(decreasing=T) %>%
   #    round()
   #
   # pred_reports$CV$prob_scaled['DO20780',] %>%
   #    sort(decreasing=T) %>%
   #    round(3)
   #
   # round(sort(pred_reports$holdout$prob_scaled['XXXXXXXX',]) * 100)
   # round(sort(pred_reports$CV$prob_scaled['DO20780',]) * 100)

   # correct_stats <- lapply(1:ncol(topn_correct), function(i_topn){
   #    #i_topn=1
   #    correct_stats.topn <- lapply(c('All',uniq_classes), function(j_class){
   #       #i_topn=1
   #       #j_class='HeadAndNeck_Other'
   #
   #       if(j_class=='All'){
   #          is_actual_class <- rep(TRUE, length(actual))
   #       } else {
   #          is_actual_class <- actual==j_class
   #       }
   #
   #       j_df <- data.frame(
   #          class=j_class,
   #          is_actual_class=is_actual_class,
   #          topn=i_topn,
   #          is_correct_pred=as.integer(topn_correct[,i_topn]),
   #          mccs=mccs
   #       )
   #       j_df$is_correct_pred <- as.integer(topn_correct[,i_topn])
   #       j_df$is_correct_pred[]
   #
   #       ## Calculate fraction correctly predicted
   #       j_df <- j_df[order(j_df$mccs, decreasing=T),]
   #       n_correct <- sum(actual_class)
   #       j_df$n_correct <- sapply(j_df$is_correct_pred, function(bool){
   #          if(!bool){ n_correct <<- n_correct - 1 }
   #          return(n_correct)
   #       })
   #       j_df$frac_correct <- j_df$n_correct/sum(actual_class)
   #
   #       # j_df <- j_df[order(j_df$mccs, decreasing=F),]
   #       # #subset(j_df, actual_class==F)
   #       # j_df$n_correct <- cumsum(j_df$is_correct_pred)
   #       # j_df$frac_correct <- j_df$n_correct/sum(actual_class)
   #
   #       return(j_df)
   #    })
   #    correct_stats.topn <- do.call(rbind, correct_stats.topn)
   #    return(correct_stats.topn)
   # })
   # correct_stats <- do.call(rbind, correct_stats)

   correct_stats <- lapply(1:ncol(topn_correct), function(i_topn){
      #i_topn=1
      correct_stats.topn <- lapply(c('All',uniq_classes), function(j_class){
         #i_topn=1
         #j_class='HeadAndNeck_Other'

         if(j_class=='All'){
            sel_samples <- rep(TRUE, length(actual))
         } else {
            sel_samples <- actual==j_class
         }

         j_df <- data.frame(
            class=j_class,
            actual_class=actual[sel_samples],
            topn=i_topn,
            is_correct_pred=topn_correct[sel_samples,i_topn],
            mccs=mccs[sel_samples]
         )

         ## Calculate fraction correctly predicted
         j_df <- j_df[order(j_df$mccs, decreasing=T),]
         n_correct <- nrow(j_df)
         j_df$n_correct <- sapply(j_df$is_correct_pred, function(bool){
            if(!bool){ n_correct <<- n_correct - 1 }
            return(n_correct)
         })
         j_df$frac_correct <- j_df$n_correct / nrow(j_df)

         ##########
         # ## Calculate fraction correctly predicted
         # j_df <- j_df[order(j_df$mccs, decreasing=F),]
         # j_df$n_correct <- cumsum(j_df$is_correct_pred)
         # j_df$frac_correct <- j_df$n_correct / nrow(j_df)
         #
         # ## Add values when MCCS is 0 or 1
         # j_df <- rbind(
         #    within(j_df[1,], { mccs <- 0; n_correct <- 0; frac_correct <- 0 }),
         #    j_df,
         #    within(j_df[nrow(j_df),], { mccs <- 1 })
         # )
         #
         # ##
         # j_df <- j_df[order(j_df$mccs, decreasing=T),]
         return(j_df)
      })
      correct_stats.topn <- do.call(rbind, correct_stats.topn)
      return(correct_stats.topn)
   })
   correct_stats <- do.call(rbind, correct_stats)

   ggplot(subset(correct_stats, topn<=3), aes(x=mccs, y=frac_correct, color=as.factor(topn))) +
      facet_wrap(~class) +
      geom_path(size=0.5) +
      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      theme_bw()


   #####
   pd_pre <- lapply(1:4, function(topn){
      #topn=1
      topn_string <- paste0('Top ',topn)

      ##
      coords <- data.frame(
         mccs,
         is_correct_pred=topn_correct[,topn]
      )
      coords <- coords[order(coords$mccs, decreasing=T),]
      coords$topn <- topn_string

      ##
      curve <- isoReg(x=coords$mccs, y=coords$is_correct_pred)
      curve$topn <- topn_string

      if(F){
         predict.isoReg(curve, 0.825)
         dim(subset(coords, mccs>=0.825 & topn=='Top 3'))
         dim(subset(coords, topn=='Top 3'))
      }

      list(coords=coords, curve=curve)
   })

   coords <- do.call(rbind, lapply(pd_pre,`[[`,'coords'))
   curves <- do.call(rbind, lapply(pd_pre,`[[`,'curve'))

   ggplot() +
      geom_point(data=coords, mapping=aes(x=mccs, y=as.integer(is_correct_pred), color=is_correct_pred), shape=21, alpha=0.5) +
      geom_density(data=coords, mapping=aes(x=mccs, y=..scaled.., color=is_correct_pred)) +
      #geom_histogram(data=coords, mapping=aes(x=mccs, color=is_correct_pred)) +
      geom_line(data=curves, mapping=aes(x=x, y=y)) +
      facet_wrap(~topn) +

      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      labs(
         x='Multiclass combined score',
         y='Probability of correct pred',
         title='Black line: combined score --> prob of correct pred calibration curve\nRed/blue lines: density of samples (dots)',
         color='Has correct prediction\nin top-n'
      ) +

      theme_bw()



   confusion <- confusionMatrixFloat(
      actual.bool=topn_correct[,1],
      probs.predicted=mccs
   )

   confusion <- as.data.frame(confusion)
   #confusion$frac_correct <- with( confusion, (tp+tn)/(tp+fp+tn+fn) )
   #confusion$frac_correct <- with( confusion, (fn+tn)/(tp+fp+tn+fn) )
   confusion$frac_correct <- with( confusion, fn/(tp+fp+tn+fn) )
   confusion <- confusion[order(confusion$frac_correct, decreasing=T),]
   #head(subset(confusion, frac_correct<0.8))

   ggplot(confusion, aes(x=prob, y=frac_correct)) +
      geom_line() +
      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      theme_bw()
}


confusionMatrixFloat <- function(actual.bool, probs.predicted){

   # if(F){
   #    ##
   #    i='Thyroid'
   #    df <- data.frame(
   #       class=i,
   #       prob=probs[,i],
   #       actual=actual==i#,
   #       #predicted=predicted==i
   #    )
   #    #head(df, 100)
   #
   #    ##
   #    i='Thyroid'
   #    actual.bool = actual==i
   #    probs.predicted = probs[,i]
   # }

   df <- data.frame(
      prob=probs.predicted,
      actual=actual.bool
   )

   pos_total <- sum(df$actual==TRUE)
   neg_total <- sum(df$actual==FALSE)

   ## actual==TRUE; TP and FN --------------------------------
   ## No. of TP = no. of TRUE samples with a prob above the current prob
   ## As prob goes from 1 to 0 (actually: highest to lowest), no. of TRUE samples with a prob above the current prob increases incrementally

   df_pos <- subset(df, actual==TRUE)
   df_pos <- df_pos[order(df_pos$prob, decreasing=T),]

   current_prob <- 1
   current_tp <- 0
   df_pos$tp <- sapply(df_pos$prob, function(prob_value){
      if(prob_value<=current_prob){
         current_tp <<- current_tp + 1
         current_prob <<- prob_value
      }
      return(current_tp)
   })

   ## FN is the inverse of TP
   df_pos$fn <- pos_total - df_pos$tp

   ##
   df_pos$not_dup <- !duplicated(df_pos$prob, fromLast=T)

   ## De-duplicate; Add 0 and 1 prob rows
   stats_pos <- subset(df_pos, not_dup, c(prob, tp, fn)); rownames(stats_pos) <- NULL

   if(stats_pos$prob[1]!=1){
      stats_pos <- rbind(
         within(stats_pos[1,],{ prob <- 1; tp <- 0; fn <- pos_total }),
         stats_pos
      )
   }

   if(tail(stats_pos$prob,1)!=0){
      stats_pos <- rbind(
         stats_pos,
         within(stats_pos[1,],{ prob <- 0; tp <- pos_total; fn <- 0 })
      )
   }

   #head(stats_pos)
   #tail(stats_pos)

   ## actual==FALSE; FP and TN --------------------------------
   ## No. of FP = no. of FALSE samples with a prob above the current prob
   ## As prob goes from 1 to 0 (actually: highest to lowest), no. of FALSE samples with a prob above the current prob increases incrementally

   df_neg <- subset(df, actual==FALSE)
   df_neg <- df_neg[order(df_neg$prob, decreasing=T),]

   ##
   current_prob <- 1
   current_fp <- 0
   df_neg$fp <- sapply(df_neg$prob, function(prob_value){
      if(prob_value <= current_prob){
         current_fp <<- current_fp + 1
         current_prob <<- prob_value
      }
      return(current_fp)
   })

   ## TN is the inverse of FP
   df_neg$tn <- neg_total - df_neg$fp

   ##
   df_neg$not_dup <- !duplicated(df_neg$prob, fromLast=T)

   #tail(df_neg, 25)
   #head(df_neg, 25)

   ## De-duplicate; Add 0 and 1 prob rows
   stats_neg <- subset(df_neg, not_dup, c(prob, fp, tn)); rownames(stats_neg) <- NULL

   if(stats_neg$prob[1]!=1){
      stats_neg <- rbind(
         within(stats_neg[1,],{ prob <- 1; fp <- 0; tn <- neg_total }),
         stats_neg
      )
   }

   if(tail(stats_neg$prob,1)!=0){
      stats_neg <- rbind(
         stats_neg,
         within(stats_neg[1,],{ prob <- 0; fp <- neg_total; tn <- 0 })
      )
   }
   #head(stats_neg)
   #tail(stats_neg)

   ## Interpolate t/f p/n for all existing probs --------------------------------
   all_probs <- sort(unique(c(stats_pos$prob, stats_neg$prob)))

   cbind(
      prob=all_probs,
      tp=approx(x=stats_pos$prob, y=stats_pos$tp, xout=all_probs)$y,
      fn=approx(x=stats_pos$prob, y=stats_pos$fn, xout=all_probs)$y,
      fp=approx(x=stats_neg$prob, y=stats_neg$fp, xout=all_probs)$y,
      tn=approx(x=stats_neg$prob, y=stats_neg$tn, xout=all_probs)$y
   )
}

