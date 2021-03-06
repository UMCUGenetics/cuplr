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
#' @param n.tasks.per.node Number of threads to use. If NULL, the number of threads will be the
#' number of classes
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
   time='2:00:00',mem='16G', n.tasks.per.node=NULL,
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

   if(is.null(n.tasks.per.node)){ n.tasks.per.node <- n_classes }

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
#SBATCH --ntasks-per-node=${n.tasks.per.node}
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
      job_script <- gsub2("${n.tasks.per.node}",n.tasks.per.node)

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
#' Gathers results from CV fold dirs and creates summary plots
#'
#' @param cv.out.dir CV dir. Output path from `spawnCvJobs()`
#' @param out.path Path to write output
#' @param pattern Name/pattern of the rds file outputted by the training function
#' @param verbose Show progress messages?
#'
#' @return Writes merged CV data and plots to out.dir
#' @export
#'
gatherCvOutput <- function(cv.out.dir, out.path, pattern='^test_set_report.rds$',verbose=T ){

   if(!file.exists(out.path)){
      if(verbose){ message('Reading CV output...') }
      report_paths <- list.files(path=cv.out.dir, pattern=pattern, recursive=T, full.names=T)
      reports <- lapply(report_paths, function(i){
         if(verbose){ message('  ',i) }
         readRDS(i)
      })

      if(verbose){ message('Merging CV output...') }
      reports_merged <- combineLists(reports, exclude.objects='imp', verbose=verbose)

      if(verbose){ message('Aggregating importances...') }
      reports_merged$imp <- aggregateMatrixList(lapply(reports,`[[`,'imp'), as.matrix=T)

      if(verbose){ message('Gathering fold numbers...') }
      fold_n_samples <- sapply(reports, function(i){
         length(i$class_pred)
      })
      reports_merged$fold_num <- rep(1:length(reports), fold_n_samples)

      if(verbose){ message('Saving merged CV results: ', out.path) }
      saveRDS(reports_merged, out.path)
   } else {
      message('Merged CV results exist: ', out.path)
   }

}






