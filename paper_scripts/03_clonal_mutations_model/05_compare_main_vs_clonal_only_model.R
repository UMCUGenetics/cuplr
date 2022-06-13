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
wd <-  paste0(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21b_DR104-update5_pcawgSoftFilt2_clonal/'))
model_dir <- wd
#model_dir <- paste0(paste0(wd,'/seed_models/01/'))
#out_dir <- paste0(model_dir,'/report/')
out_dir <- paste0(wd,'/report/')
dir.create(out_dir, showWarnings=F)

#model <- cacheAndReadData(paste0(model_dir,'/final/model.rds'), overwrite=T)
#features <- cacheAndReadData(paste0(model_dir,'/features/features_allSamples.rds'), overwrite=T)
metadata <- read.delim(paste0(wd,'/features/metadata_training.txt'))


## Pred reports --------------------------------
## CV set
gatherCvOutput(
   cv.out.dir=paste0(model_dir,'/cv_out/'),
   out.path=paste0(model_dir,'/test_set_report.rds')
)
pred_reports.cv <- cacheAndReadData(paste0(model_dir,'/test_set_report.rds'), overwrite=T)

pred_reports <- list(CV=pred_reports.cv)

pred_reports <- lapply(pred_reports, function(i){
   i$class_actual <- as.factor(
      metadata$cancer_type[ match(rownames(i$prob), metadata$sample_id) ]
   )
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
   df <- summary.predReport(pred_reports[[i]], top.n.classes=3, top.n.feat=5, top.n.feat.classes=3, prob.type='prob_scaled')
   data.frame(group=i, df)
}))
#View(subset(pred_summ, group=='CUP'))

if(WRITE_OUTPUT){
   write.table(pred_summ, paste0(out_dir,'/pred_summ.txt'), sep='\t', quote=F, row.names=F)
}

## Compare performance with main model ================================
## Load data --------------------------------
l_pred_summ <- list()

l_pred_summ$all_muts <- (function(){
   df <- read.delim(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/report/pred_summ.txt'))
   subset(df, group=='CV')
})()

l_pred_summ$clonal_muts <- subset(pred_summ, group=='CV')

## Plot data --------------------------------
pd <- lapply(names(l_pred_summ), function(i){
   df <- l_pred_summ[[i]]
   out <- calcFracCorrect(
      as.factor(df$actual_class),
      as.factor(df$pred_class.1)
   )
   out$group <- i
   return(out)
})
pd <- do.call(rbind, pd)

pd$value_string <- with(pd, paste0(
   round(frac_correct,2),' (',n_correct,'/',n_total,')'
))

## Plot --------------------------------
p <- ggplot(pd, aes(x=group, y=frac_correct)) +
   facet_wrap(~class) +
   #geom_bar(aes(fill=group), stat='identity', position=position_dodge(width=0.9)) +
   geom_bar(aes(fill=group), stat='identity', color='black', size=0.3) +
   geom_text(aes(label=value_string, y=0.04), angle=90, hjust=0, vjust=0.5, size=2.7) +
   ylab('Fraction of samples correctly predicted') +
   scale_fill_discrete(name='Which mutations') +
   theme_bw() +
   theme(
      #axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
   )

pdf(paste0(out_dir,'/all_vs_clonal_cuplr.pdf'), 11, 8)
plot(p)
dev.off()
