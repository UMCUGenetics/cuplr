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
   k=20, seed=1, rm.path.prefix='/Users/lnguyen',
   time='2:00:00',mem='16G',
   verbose=T
){
   if(F){
      #train.data.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.06d_probWeighRf_balanceClasses_slurm/features/features.rds'
      train.data.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10c_rmdTad_DR104update/features/features.txt.gz'
      train.script.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10c_rmdTad_DR104update/do_train.R'
      colname.response='response'
      k=20
      seed=1
      cv.out.dir=paste0(out_dir,'/cv_out/')
      rm.path.prefix='/Users/lnguyen'
   }

   set.seed(seed)

   if(verbose){ message('Creating jobs @: ', cv.out.dir) }
   dir.create(cv.out.dir, showWarnings=F, recursive=T)

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
Rscript ${train.script.path} ${train.data.path} ${fold_indexes_path} ${out_path} ${seed} && touch ${done_path}
EOF
else
echo Skipping ${fold_name}. Done file exists: ${done_path}
fi
")
      job_script <- gsub("${fold_name}",fold_name, job_script, fixed=T)
      job_script <- gsub("${fold_dir}",fold_dir, job_script, fixed=T)

      job_script <- gsub("${time}",time, job_script, fixed=T)
      job_script <- gsub("${mem}",mem, job_script, fixed=T)

      job_script <- gsub("${n_classes}",n_classes, job_script, fixed=T)

      job_script <- gsub("${train.script.path}",train.script.path, job_script, fixed=T)
      job_script <- gsub("${train.data.path}",train.data.path, job_script, fixed=T)
      job_script <- gsub("${fold_indexes_path}",fold_indexes_path, job_script, fixed=T)
      job_script <- gsub("${out_path}",out_path, job_script, fixed=T)
      job_script <- gsub("${seed}",seed, job_script, fixed=T)
      job_script <- gsub("${done_path}",done_path, job_script, fixed=T)

      job_script <- gsub(rm.path.prefix,'', job_script, fixed=T)

      job_script_path <- paste0(fold_dir,'/job.sh')
      file_conn <- file(job_script_path)
      writeLines(job_script, file_conn)
      close(file_conn)
   }
}

####################################################################################################
#' Gathers results from CV fold dirs and creates summary plots
#'
#' @param cv.out.dir CV dir. Output path from `spawnCvJobs()`
#' @param out.dir Output dir for this function
#' @param pattern Name/pattern of the rds file outputted by the training function
#' @param verbose Show progress messages?
#'
#' @return Writes merged CV data and plots to out.dir
#' @export
#'
gatherCvOutput <- function(
   cv.out.dir, out.dir=paste0(cv.out.dir,'/../'),
   pattern='model.rds', verbose=T
){
   #cv.out.dir='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.06d_probWeighRf_balanceClasses_slurm/cv_out/'
   #cv.out.dir='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/training/models/0.06e_balanceClasses/cv_out/'
   output_paths <- list.files(path=cv.out.dir, pattern=pattern, recursive=T, full.names=T)

   out_files <- paste0(
      out.dir,
      c('/test_set.rds','/imp.rds')
   )
   names(out_files) <- c('test_set','imp')


   if(!all(file.exists(out_files))){
      if(verbose){ message('Reading output from CV folds...') }

      l_test_set <- list()
      l_imp <- list()

      for(i in 1:length(output_paths)){
         #i=1
         if(verbose){ message('> ',output_paths[[i]]) }
         output <- readRDS(output_paths[[i]])
         l_test_set[[i]] <- output$test_set
         l_imp[[i]] <- output$imp
      }

      test_set <- list(
         probabilities=do.call(rbind, lapply(l_test_set,`[[`,'probabilities')),
         predicted=unlist(lapply(l_test_set,`[[`,'predicted')),
         actual=unlist(lapply(l_test_set,`[[`,'actual'))
      )

      imp <- aggregateMatrixList(l_imp, as.matrix=T)

      saveRDS(test_set, out_files['test_set'])
      saveRDS(imp, out_files['imp'])
   } else {
      if(verbose){ message('Reading merged CV data...') }
      test_set <- readRDS(out_files['test_set'])
      imp <- readRDS(out_files['imp'])
   }

   plots_dir <- paste0(out.dir,'/plots/')
   dir.create(plots_dir, showWarnings=F)

   if(verbose){ message('Plotting imp barplots...') }
   pdf(paste0(plots_dir,'/imp_barplots.pdf'), 16, 10)
   plot( plotTopFeatures(imp, top.n=40, infer.feature.type=T, n.col=4, feature.type.colors=NULL) )
   dev.off()

   # if(verbose){ message('Plotting feat imp heatmap...') }
   # pdf(paste0(plots_dir,'/imp_heatmap.pdf'), 17, 8)
   # plot( plotFeatureImpHeatmap(imp, top.n=150) )
   # dev.off()

   if(verbose){ message('Plotting perf heatmap...') }
   pdf(paste0(plots_dir,'/perf_heatmap.pdf'), 12, 10)
   suppressWarnings({
      plot(plotPerfHeatmap(
         test_set$actual, test_set$predicted, show.weighted.mean=T,
         rel.heights=c(0.3, 0.12, 1)
      ))
   })
   dev.off()

   # if(verbose){ message('Plotting false negative rates...') }
   # pdf(paste0(plots_dir,'/fnr_heatmap.pdf'), 11, 8)
   # plot( plotFnr(test_set$actual, test_set$prob) )
   # dev.off()

}






