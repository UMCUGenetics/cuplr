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
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor'))
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/statsExtra/'))
devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/'))

## Load data ================================
## Gather CV predictions --------------------------------
model_dir <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21c_featGroupExcl/')
#model_dir <- paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21c_featGroupExcl/seed_models/17/')
pred_reports_dir <- paste0(model_dir,'/pred_reports/')
dir.create(pred_reports_dir, showWarnings=F)

cv_out_subdirs <- list.dirs(paste0(model_dir,'/cv_out/'), recursive=F)
cv_out_subdirs <- cv_out_subdirs[ grepl('^\\d+', basename(cv_out_subdirs)) ]

for(i in 1:length(cv_out_subdirs)){
   #i=1
   gatherCvOutput(
      cv.out.dir=cv_out_subdirs[[i]],
      out.path=paste0(pred_reports_dir,'/',basename(cv_out_subdirs[[i]]),'.test_set_report.rds')
   )
}

## Load CV predictions --------------------------------
pred_reports_paths <- list.files(pred_reports_dir, full.names=T)
names(pred_reports_paths) <- basename(pred_reports_paths)
names(pred_reports_paths) <- sub('[.].+','',names(pred_reports_paths))
names(pred_reports_paths) <- sub('\\d+_','',names(pred_reports_paths))

pred_reports_paths <- c(
   none=paste0(base_dir,'/CUPs_classifier/processed/cuplr/cuplr/models/0.21a_DR104-update5_pcawgSoftFilt2/seed_models/01/test_set_report.rds'),
   pred_reports_paths
)

pred_reports <- lapply(pred_reports_paths, function(i){
   report <- cacheAndReadData(i, overwrite=T)
   class(report) <- c(class(report),'predReport')
   return(report)
})

## Prob calibration ================================
## Isotonic regression
cal_curves <- lapply(pred_reports, function(i){
   probCalCurves(
      actual=i$class_actual,
      probs=i$prob,
      method='isotonic'
   )
})

## Apply calibration curve
pred_reports <- lapply(names(pred_reports), function(i){
   report <- pred_reports[[i]]
   report$prob_scaled <- probCal(report$prob, cal.curves=cal_curves[[i]])

   report$class_pred <- colnames(report$prob_scaled)[ max.col(report$prob_scaled) ]
   report$class_pred <- factor(report$class_pred, colnames(report$prob_scaled))

   return(report)
})
names(pred_reports) <- names(pred_reports_paths)

## Perf ================================
if(WRITE_OUTPUT){

   feature_type_colors <- c(
      #none="#FB8072",## red
      none="white",
      FEATURE_TYPE_COLORS
   )

   perf <- (function(){
      df <- lapply(names(pred_reports), function(i){
         #i='none'
         report <- pred_reports[[i]]

         out <- calcFracCorrect(actual=report$class_actual, predicted=report$class_pred)
         out <- cbind(excl_feat_group=i, out)
         return(out)
      })
      df <- do.call(rbind, df)
      df$n_incorrect <- df$n_total - df$n_correct

      ## Assign factors
      df$excl_feat_group <- factor(df$excl_feat_group, unique(df$excl_feat_group))
      df$class <- factor(df$class, unique(df$class))


      df$excl_feat_group <- factor( df$excl_feat_group, names(feature_type_colors) )
      df$excl_feat_group <- droplevels(df$excl_feat_group)

      ## Labels
      df$facet_label <- paste0(df$class,' (',df$n_total,')')
      df$bar_label <- paste0(round(df$frac_correct,3),' (',df$n_correct,')')

      return(df)
   })()

   ##
   p_perf <- ggplot(perf, aes(x=excl_feat_group, y=frac_correct)) +
      facet_wrap(.~facet_label, scales='free_y', ncol=4) +

      ## Main bars
      geom_bar(aes(fill=excl_feat_group), stat='identity', color='black', size=0.3) +

      ## Line marking model with no feature types removed
      geom_hline(
         data=subset(perf, excl_feat_group=='none'), mapping=aes(yintercept=frac_correct),
         size=0.3, color='red'
      ) +

      ## Bar labels
      geom_text(aes(label=bar_label, y=0.02), angle=90, hjust=0, vjust=0.5, size=2.5) +

      scale_fill_manual(values=feature_type_colors, name='Excluded\nfeature type') +
      xlab('Excluded feature type') +
      ylab('Recall') +

      theme_bw() +
      theme(
         panel.grid.minor=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
      )

   pdf(paste0(model_dir,'/acc_feat_excl.pdf'), 11, 11)
   plot(p_perf)
   dev.off()
}
