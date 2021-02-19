#' Spawn CV jobs for submission with SLURM
#'
#' @param df A dataframe of features and the response variable
#' @param train.data.path Path to rds file. A dataframe of features and the response variable
#' @param train.script.path Path to training R script. The following input arguments are required in
#' within the script:
#' args <- commandArgs(trailingOnly=T)
#' training_data <- args[1]
#' fold_indexes <- args[2]
#' out_path <- args[3]
#' seed <- args[4]
#' @param cv.out.dir CV output dir. Jobs, input data and output data will be stored here
#' @param colname.response The column name of the response variable (i.e. training labels).
#' @param k Number of cross validation folds.
#' @param seed Random seed as in integer.
#' @param rm.path.prefix Prefix to remove from all paths in the job script
#' @param trailing.args A string specifying trailing args in Rscript command:
#' Rscript ${train.script.path} ${train.data.path} ${fold_indexes_path} ${out_path} ${seed} ${trailing.args}
#' @param time A string specifying SLURM run time
#' @param mem A string specifying SLURM memory
#' @param verbose Show progress messages?
#'
#' @return Writes jobs to cv.out.dir
#' @export
#'
spawnCvJobs <- function(
   train.data.path=NULL, df=NULL,
   train.script.path,
   cv.out.dir=paste0(dirname(train.script.path),'/cv_out/'),
   colname.response='response',
   k=20, seed=1,
   rm.path.prefix='/Users/lnguyen', trailing.args='',
   time='2:00:00',mem='16G',
   verbose=T
){
   if(F){
      #train.data.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.06d_probWeighRf_balanceClasses_slurm/features/features.rds'
      wd='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.13a_HMF_PCAWG/'
      train.data.path=paste0(wd,'/features/training_labels.txt')
      train.script.path=paste0(wd,'/do_train.R')
      colname.response='response'
      k=15
      seed=1
      cv.out.dir=paste0(wd,'/cv_feature_excl/')
      rm.path.prefix='/Users/lnguyen'
      verbose=T
   }

   set.seed(seed)

   if(verbose){ message('Creating jobs @: ', cv.out.dir) }
   dir.create(cv.out.dir, showWarnings=F, recursive=T)

   ## --------------------------------
   if(!is.data.frame(df)){
      if(verbose){ message('Reading training data: ', train.data.path) }
      if(grepl('.rds$',train.data.path)){
         df <- readRDS(train.data.path)
      } else {
         df <- read.delim(train.data.path, check.names=F)
      }
   }

   n_classes <- length(unique(df[,colname.response]))

   if(verbose){ message('Getting fold indexes...') }
   folds <- mltoolkit::createCvTrainTestSets(df, stratify.by.col=colname.response, k=k, return.data=F)

   ## Helper functions --------------------------------
   writeString <- function(string,path){
      file_conn <- file(path)
      writeLines(string, file_conn)
      close(file_conn)
   }

   ## --------------------------------
   if(verbose){ message('Writing jobs...') }
   for(i in 1:length(folds)){
      #i=1
      fold_name <- formatC(i,width=nchar(length(folds)),format="d",flag="0")
      if(verbose){ message('> Fold: ',fold_name) }

      fold_dir <- paste0(cv.out.dir,'/',fold_name,'/')
      dir.create(fold_dir, showWarnings=F, recursive=F)

      fold_indexes_path <- paste0(fold_dir,'/fold_indexes.rds')
      saveRDS(folds[[i]], fold_indexes_path)

      out_path <- paste0(fold_dir,'/model.rds')
      done_path <- paste0(fold_dir,'/job.done')

      job_script <- paste0("#!/bin/bash
#SBATCH --job-name=cv_${fold_name}
#SBATCH --output=${fold_dir}/slurm.out
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH --ntasks-per-node=${n_classes}
if [[ ! -f ${done_path} ]]; then
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript ${train.script.path} ${train.data.path} ${fold_indexes_path} ${out_path} ${seed} ${trailing.args} && touch ${done_path}
EOF
else
echo Skipping ${fold_name}. Done file exists: ${done_path}
fi
")
      gsub2 <- function(var.string, replace.string){
         gsub(var.string, replace.string, job_script, fixed=T)
      }

      ## SBATCH header
      job_script <- gsub2("${fold_name}",fold_name)
      job_script <- gsub2("${fold_dir}",fold_dir)
      job_script <- gsub2("${time}",time)
      job_script <- gsub2("${mem}",mem)
      job_script <- gsub2("${n_classes}",n_classes)

      ## Main
      job_script <- gsub2("${train.script.path}",train.script.path)
      job_script <- gsub2("${train.data.path}",train.data.path)
      job_script <- gsub2("${fold_indexes_path}",fold_indexes_path)
      job_script <- gsub2("${out_path}",out_path)
      job_script <- gsub2("${seed}",seed)
      job_script <- gsub2("${done_path}",done_path)
      job_script <- gsub2("${trailing.args}",trailing.args)

      job_script <- gsub(rm.path.prefix,'', job_script)

      writeString( job_script, paste0(fold_dir,'/job.sh') )
   }

   ## --------------------------------
   submit_script <- paste0("#!/bin/bash
for i in ${cv.out.dir}/*/; do
	if [[ ! -f $i/job.done ]]; then
		sbatch $i/job.sh
	fi
done
")

   submit_script <- gsub('${cv.out.dir}', cv.out.dir, submit_script, fixed=T)
   submit_script <- gsub(rm.path.prefix,'', submit_script)

   writeString( submit_script, paste0(cv.out.dir,'/submit_jobs.sh') )
}

####################################################################################################
#' Combine results from multiple prediction reports into one list object
#'
#' @param reports A list of reports. Each report (a list) must have the names: probs_raw,
#' probs_adjusted, responses_pred, responses_actual, feat_contrib, imp
#' @param verbose Show progress messages?
#'
#' @return A list
#' @export
#'
mergePredReports <- function(reports, verbose=T){

   reports_merged <- list()

   if(verbose){ 'Merging probs_raw...' }
   reports_merged$probs_raw <- do.call(rbind, lapply(reports,`[[`,'probs_raw'))

   if(verbose){ 'Merging probs_adjusted...' }
   reports_merged$probs_adjusted <- do.call(rbind, lapply(reports,`[[`,'probs_adjusted'))

   if(verbose){ 'Merging predicted responses...' }
   reports_merged$responses_pred <- structure(
      unlist(lapply(reports,`[[`,'responses_pred')),
      names=unlist(lapply(reports,function(i){ names(i$responses_pred) }))
   )

   if(verbose){ 'Merging actual responses...' }
   reports_merged$responses_actual <- structure(
      unlist(lapply(reports,`[[`,'responses_actual')),
      names=unlist(lapply(reports,function(i){ names(i$responses_actual) }))
   )

   if(verbose){ 'Merging feature contributions...' }
   reports_merged$feat_contrib <- do.call(rbind, lapply(reports,`[[`,'feat_contrib'))
   reports_merged$feat_contrib <- subset(reports_merged$feat_contrib, contrib>0)

   if(verbose){ 'Merging feature importances...' }
   reports_merged$imp <- aggregateMatrixList(lapply(reports,`[[`,'imp'), as.matrix=T)

   return(reports_merged)
}

####################################################################################################
#' Gathers results from CV fold dirs and creates summary plots
#'
#' @param cv.out.dir CV dir. Output path from `spawnCvJobs()`
#' @param out.dir Output dir for this function
#' @param pattern Name/pattern of the rds file outputted by the training function
#' @param mk.plots Make performance and feature importance plots?
#' @param verbose Show progress messages?
#'
#' @return Writes merged CV data and plots to out.dir
#' @export
#'
gatherCvOutput <- function(
   cv.out.dir,
   out.dir=paste0(cv.out.dir,'/../'),
   pattern='^test_set_report.rds$',
   mk.plots=T,
   verbose=T
){
   if(F){
      cv.out.dir='/Users/lnguyen//hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.13a_HMF_PCAWG/cv_out/'
      out.dir=paste0(cv.out.dir,'/../')
      pattern='^test_set_report.rds$'
      verbose=T
   }

   if(verbose){ message('Reading CV output...') }
   report_paths <- list.files(path=cv.out.dir, pattern=pattern, recursive=T, full.names=T)
   reports <- lapply(report_paths, function(i){
      if(verbose){ message('  ',i) }
      readRDS(i)
   })

   ## Combine results from cv folds --------------------------------
   reports_merged_path <- paste0(out.dir,'/test_set_report.rds')

   if(!file.exists(reports_merged_path)){
      if(verbose){ message('Merging CV results...') }
      reports_merged <- mergePredReports(reports)

      fold_n_samples <- sapply(reports, function(i){
         length(i$responses_pred)
      })
      reports_merged$fold_num <- rep(1:length(reports), fold_n_samples)

      #rep(c(1,2),c(2,3))

      saveRDS(reports_merged, reports_merged_path)
   } else {
      if(verbose){ message('Loading merged CV results: ', reports_merged_path) }
      reports_merged <- readRDS(reports_merged_path)
   }

   ## Plots --------------------------------
   if(mk.plots){
      plots_dir <- paste0(out.dir,'/plots/')
      dir.create(plots_dir, showWarnings=F)

      if(verbose){ message('Plotting perf heatmap...') }
      pdf(paste0(plots_dir,'/perf_heatmap.pdf'), 12, 10)
      suppressWarnings({
         plot(plotPerfHeatmap(
            reports_merged$responses_actual, reports_merged$responses_pred, show.weighted.mean=T,
            rel.heights=c(perf=0.3, counts=0.15, heatmap=1)
         ))
      })
      dev.off()

      if(verbose){ message('Plotting imp barplots...') }
      pdf(paste0(plots_dir,'/imp_barplots.pdf'), 16, 10)
      plot( plotTopFeatures(reports_merged$imp, top.n=40, infer.feature.type=T, n.col=4, feature.type.colors=NULL) )
      dev.off()
   }
}






