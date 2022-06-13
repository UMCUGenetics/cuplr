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
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/statsExtra/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))

## Load CV data ================================
## Gather CV predictions --------------------------------
model_dir <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21d_HMF_or_PCAWG_train/')
pred_reports_dir <- paste0(model_dir,'/pred_reports/')
dir.create(pred_reports_dir, showWarnings=F)

cv_out_subdirs <- list.dirs(paste0(model_dir,'/cv_out/'), recursive=F)
for(i in 1:length(cv_out_subdirs)){
   #i=1
   gatherCvOutput(
      cv.out.dir=cv_out_subdirs[[i]],
      out.path=paste0(pred_reports_dir,'/',basename(cv_out_subdirs[[i]]),'.test_set_report.rds')
   )
}

## Reload CV predictions --------------------------------
pred_reports_paths <- list.files(pred_reports_dir, full.names=T)
names(pred_reports_paths) <- basename(pred_reports_paths)
names(pred_reports_paths) <- sub('[.].+','',names(pred_reports_paths))
names(pred_reports_paths) <- paste0(names(pred_reports_paths),'_trained')

pred_reports.cv <- lapply(pred_reports_paths, function(i){
   report <- cacheAndReadData(i, overwrite=T)
   class(report) <- c(class(report),'predReport')
   return(report)
})

## Prob calibration ================================
## Isotonic regression
cal_curves <- lapply(pred_reports.cv, function(i){
   probCalCurves(
      actual=i$class_actual,
      probs=i$prob,
      method='isotonic'
   )
})

# ## Apply calibration curve
# pred_reports.cv <- lapply(names(pred_reports.cv), function(i){
#    report <- pred_reports.cv[[i]]
#    report$prob_scaled <- probCal(report$prob, cal.curves=cal_curves[[i]])
#
#    report$class_pred <- colnames(report$prob_scaled)[ max.col(report$prob_scaled) ]
#    report$class_pred <- factor(report$class_pred, colnames(report$prob_scaled))
#
#    return(report)
# })
# names(pred_reports.cv) <- names(pred_reports_paths)

## Predict ================================
## Load data --------------------------------
## Models
final_model_paths <- list.dirs(paste0(model_dir,'/final/'), recursive=F)
final_model_paths <- paste0(final_model_paths,'/model.rds')

final_models <- lapply(final_model_paths, function(i){ cacheAndReadData(i, overwrite=T) })
names(final_models) <- paste0(basename(dirname(final_model_paths)),'_trained')

## Misc
features <- cacheAndReadData(paste0(model_dir,'/features/features_allSamples.txt.gz'), overwrite=T)
metadata_raw <- read.delim(paste0(model_dir,'/features/metadata_training.txt'))
training_samples <- read.delim(paste0(model_dir,'/features/training_sample_selection.txt'))

metadata <- metadata_raw
metadata$is_training_sample.HMF <- training_samples$is_selected.HMF[match(metadata$sample_id, training_samples$sample)]
metadata$is_training_sample.HMF[is.na(metadata$is_training_sample.HMF)] <- FALSE
metadata$is_training_sample.PCAWG <- training_samples$is_selected.PCAWG[match(metadata$sample_id, training_samples$sample)]
metadata$is_training_sample.PCAWG[is.na(metadata$is_training_sample.PCAWG)] <- FALSE

## Predict --------------------------------
val_samples <- list(
   HMF_trained=subset(metadata, (is_training_sample | is_holdout_sample) & cohort=='PCAWG', sample_id, drop=T),
   PCAWG_trained=subset(metadata, (is_training_sample | is_holdout_sample) & cohort=='HMF', sample_id, drop=T)
)

## Main
pred_reports <- lapply(names(val_samples), function(i){
   out <- predict(
      object=final_models[[i]],
      newdata=features[val_samples[[i]],],
      calc.feat.contrib=T, type='report', top.n.pred.classes=NULL
   )
   out$class_actual <- metadata$cancer_type[ match(rownames(out$prob), metadata$sample_id) ]

   return(out)
})
names(pred_reports) <- names(val_samples)

## Prob calibration
pred_reports <- lapply(names(pred_reports), function(i){
   #i='HMF_trained'
   report <- pred_reports[[i]]

   report$prob_scaled <- probCal(report$prob, cal.curves=cal_curves[[i]])
   report$class_pred <- colnames(report$prob_scaled)[ max.col(report$prob_scaled) ]
   report$class_pred <- factor(report$class_pred, colnames(report$prob_scaled))

   return(report)
})
names(pred_reports) <- names(val_samples)

## Fix the class of the objects
pred_reports <- lapply(pred_reports, function(i){
   class(i) <- c('predReport',class(i))
   return(i)
})

## Prediction summary  --------------------------------
exist_classes <- lapply(final_models, function(i){ names(i$ensemble) })
pred_summs <- lapply(names(pred_reports), function(i){
   #i='HMF_trained'
   df <- summary(pred_reports[[i]], top.n.classes=3, top.n.feat=5, top.n.feat.classes=1, prob.type='prob_scaled')

   sel_classes <- exist_classes[[i]]
   df <- df[df$actual_class %in% sel_classes,]
   df$actual_class <- factor(df$actual_class, sel_classes)
   df$pred_class.1 <- factor(df$pred_class.1, sel_classes)

   #df$actual_class <-
   df <- data.frame(group=i, df)
   return(df)
})
names(pred_summs) <- names(pred_reports)

##
training_ct_counts <- lapply(unique(training_samples$cohort), function(i){
   df <- subset(training_samples, cohort==i & (is_selected.HMF | is_selected.PCAWG))
   ct_counts <- unclass(table(df$cancer_type))
   ct_counts <- c(All=sum(ct_counts), ct_counts)
   return(ct_counts)
})
names(training_ct_counts) <- paste0(unique(training_samples$cohort),'_trained')

acc <- lapply(names(pred_summs), function(i){
   #='HMF_trained'
   pred_summ <- pred_summs[[i]]
   df <- calcFracCorrect(actual=pred_summ$actual_class, predicted=pred_summ$pred_class.1)
   df <- df[df$n_total!=0,] ## Remove classes with no test samples
   df$n_training <- training_ct_counts[[i]][ as.character(df$class) ]
   df$group <- i
   return(df)
})
#names(acc) <- names(pred_summs)
acc <- do.call(rbind, acc)

out_dir <- paste0(model_dir,'/report/')
dir.create(out_dir, showWarnings=F)
if(WRITE_OUTPUT){
   write.table(acc, paste0(out_dir,'/frac_correct.txt'), sep='\t', quote=F, row.names=F)
}


## Plot --------------------------------
p_acc <- (function(){
   pd <- acc

   ## Fill missing rows
   rownames(pd) <- paste0(pd$group,'__',pd$class)
   uniq_groups <- sort(unique(pd$group))
   uniq_classes <- sort(unique(pd$class))
   required_rownames <- unlist(lapply(uniq_groups, function(i){ paste0(i,'__',uniq_classes) }))
   pd <- pd[required_rownames,]
   pd$group <- sub('__.+$','',required_rownames)
   pd$class <- sub('^.+__','',required_rownames)
   rownames(pd) <- NULL

   ## Labels
   pd$label.frac_correct <- paste0( round(pd$frac_correct,2), ' (', pd$n_correct, '/', pd$n_total, ')')
   pd$label.frac_correct[is.na(pd$n_total)] <- ''
   pd$n_training[is.na(pd$n_total)] <- 0

   ## Plot
   bar_width <- 0.7
   vlines <- (2:length(unique(pd$class)))-0.5
   fill_colors <- c("#82B1D1", "#F3C8A3")
   fill_labels <- c(HMF_trained='Train on Hartwig, test on PCAWG', PCAWG_trained='Train on PCAWG, test on Hartwig')

   p.frac_correct <- ggplot(pd, aes(x=class, y=frac_correct, fill=group)) +
      geom_bar(stat='identity', position=position_dodge(width=bar_width), width=bar_width, color='black', size=0.3) +
      geom_text(aes(label=label.frac_correct, y=0.01), position=position_dodge(width=bar_width), angle=90, hjust=0, vjust=0.5, size=2.7) +
      geom_vline(xintercept=vlines, size=0.5, color='lightgrey') +
      scale_fill_manual(values=fill_colors, labels=fill_labels) +
      ylab('Recall') +
      scale_x_discrete(position='top') +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         legend.position='none',
         axis.title.x=element_blank()
      )

   p.n_training <- ggplot(pd, aes(x=class, y=n_training, fill=group)) +
      coord_cartesian(ylim=c(0,800)) +
      geom_bar(stat='identity', position=position_dodge(width=bar_width), width=bar_width, color='black', size=0.3) +
      geom_text(aes(label=n_training, y=10), position=position_dodge(width=bar_width), angle=90, hjust=0, vjust=0.5, size=2.7) +
      geom_vline(xintercept=vlines, size=0.5, color='lightgrey') +
      scale_fill_manual(values=fill_colors, labels=fill_labels) +
      ylab('Number of\ntraining samples') +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         legend.position='bottom',
         legend.direction='vertical',
         legend.title=element_blank(),
         axis.title.x=element_blank()
      )

   cowplot::plot_grid(p.frac_correct, p.n_training, axis='tblr', align='v', ncol=1, rel_heights=c(1, 1.15))
})()

pdf(paste0(out_dir,'/frac_correct.pdf'), 10, 7)
plot(p_acc)
dev.off()

