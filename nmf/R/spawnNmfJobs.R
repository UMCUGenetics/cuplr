#' Non-negative matrix factorization
#'
#' @description Spawn SLURM jobs to run NNLM::nnmf() on the HPC
#'
#' @param A A numeric matrix (rows=samples, columns=features)
#' @param k.range Integer vector specifying the range of ranks
#' @param repeats Number of repeats per rank
#' @param impute.prop Fraction of values in `A` to set to NA for calculating performance
#' @param seed Random seed
#' @param verbose Show messages? Can be 0,1,2
#'
#' @return A list containing the factorized matrices and MSE values from rank search
#' @export
#'
spawnNmfJobs <- function(
   A, parent.dir,
   k.range=1:10, repeats=10, impute.prop=0.1, seed=1, verbose=T,
   func.scripts.path=MAIN_NMF_FUNCTIONS_PATH,
   time='0:45:00', mem='2G', rm.path.prefix='/Users/lnguyen', overwrite.jobs=T
){
   # if(F){
   #    A <- contexts$snv[
   #       metadata$sample_id[metadata$cancer_type=='Breast']
   #    ,]
   #
   #    k.range=2:3
   #    repeats=3
   #    impute.prop=0.1
   #    func.scripts.path='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/commonUtils/R/nmf.main.R'
   #    parent.dir='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/signatures/04_de_novo/nmf_breast/'
   #    rm.path.prefix='/Users/lnguyen'
   #    time='0:45:00'
   #    mem='2G'
   # }

   ## Load inputs --------------------------------
   if(!dir.exists(parent.dir)){
      stop('`parent.dir` doesnt exist: ',parent.dir)
   }

   if(is.matrix(A) | is.data.frame(A)){
      if(verbose){ message('Writing input matrix to disk...') }
      input_path <- paste0(parent.dir,'/m_input.txt.gz')
      write.table(A, gzfile(input_path), sep='\t', quote=F)
   } else if(is.character(A)){
      input_path <- A
      if(!file.exists(input_path)){ stop('File does not exist: ', input_path) }
   } else {
      stop('`A` must be a matrix, dataframe, or path to a txt file')
   }

   ## Helper functions --------------------------------
   writeString <- function(string, path){
      file_conn <- file(path)
      writeLines(string, file_conn)
      close(file_conn)
   }

   ## --------------------------------
   if(verbose){ message('Writing Rscript...') }
   run_nmf_script <- sprintf("options(stringsAsFactors=F)

require(NNLM)
args <- commandArgs(trailingOnly=T)

message('## Loading functions...')
source(args[1])

iter_out_dir <- args[2]

message('## Reading input matrix...')
A <- read.delim(args[3], check.names=F)

k <- as.integer(args[4])
repeat.num <- as.integer(args[5])

impute.prop <- %s
seed <- %s

message('## Running NMF to determine performance...')
nmf_out.path <- paste0(iter_out_dir,'/nmf_out.rds')
if(!file.exists(nmf_out.path)){
   nmf_out <- runNmf(
      A=A,
      k=k,
      repeat.num=repeat.num,
      seed=seed,
      max.samples=NULL,
      impute.prop=impute.prop,
      perf.metrics=c('mse_imputed','mse_perm'),
      return.perf=F,
      verbose=2
   )
   saveRDS(nmf_out, nmf_out.path)
} else {
   nmf_out <- readRDS(nmf_out.path)
}

write.table(
   nmf_out$perf,
   paste0(iter_out_dir,'/perf.txt'),
   sep='\\t', quote=F, row.names=F
)

write.table(
   t(nmf_out$clusters),
   gzfile(paste0(iter_out_dir,'/clusters.txt.gz')),
   sep='\\t', quote=F, row.names=F, col.names=F
)

", impute.prop, seed)

   run_nmf_script_path <- paste0(parent.dir,'/run_nmf.R')
   writeString(run_nmf_script, run_nmf_script_path)

   ## Make jobs --------------------------------
   output_dir <- paste0(parent.dir,'/output/')
   dir.create(output_dir, showWarnings=F)

   if(!overwrite.jobs){
      if(verbose){ message('`overwrite.jobs` is FALSE. Skipping writing jobs') }
      return()
   }

   ## Init
   param_grid <- data.frame(
      k=rep(k.range, each=repeats),
      repeat_num=rep(1:repeats, length(k.range))
   )

   param_grid <- as.data.frame(lapply(param_grid, function(i){
      formatC(i,width=nchar(max(i)),format="d",flag="0")
   }))
   param_grid$time <- time
   param_grid$mem <- mem

   ## Make jobs for every param set
   if(verbose>=1){ message('Writing jobs') }
   for(i in 1:nrow(param_grid)){
      #i=1
      k <- param_grid$k[i]
      repeat_num <- param_grid$repeat_num[i]
      i_time <- param_grid$time[i]
      i_mem <- param_grid$mem[i]

      job_name <- paste0('k_',k,'.rep_',repeat_num,'.job')

      if(verbose>=2){ message('> ', job_name) }

      i_out_dir <- paste0(output_dir,'/',sub('[.]job$','',job_name),'/')
      dir.create(i_out_dir, showWarnings=F)

      job_script <- "#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${i_out_dir}/${job_name}.o
#SBATCH --time=${i_time}
#SBATCH --mem=${i_mem}
if [[ ! -f ${i_out_dir}/${job_name}.done ]]; then
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript ${run_nmf_script_path} ${func.scripts.path} ${i_out_dir} ${input_path} ${k} ${repeat_num} && touch ${i_out_dir}/${job_name}.done
EOF
else
echo Skipping ${job_name}
fi
"
      gsub2 <- function(var.string, replace.string){
         gsub(var.string, replace.string, job_script, fixed=T)
      }

      job_script <- gsub2("${job_name}",job_name)
      job_script <- gsub2("${i_time}",i_time)
      job_script <- gsub2("${i_mem}",i_mem)

      job_script <- gsub2("${run_nmf_script_path}",run_nmf_script_path)
      job_script <- gsub2("${func.scripts.path}",func.scripts.path)
      job_script <- gsub2("${i_out_dir}",i_out_dir)
      job_script <- gsub2("${input_path}",input_path)
      job_script <- gsub2("${k}",k)
      job_script <- gsub2("${repeat_num}",repeat_num)

      job_script <- gsub(rm.path.prefix,'', job_script)

      writeString(job_script, paste0(i_out_dir,'/',job_name))
   }

   ## Make submit script --------------------------------
   submit_script <- paste0("#!/bin/bash
for i in ${output_dir}/*/*.job; do
	if [[ ! -f ${i}.done ]]; then
		sbatch $i
	fi
done
")
   submit_script <- gsub('${output_dir}', output_dir, submit_script, fixed=T)
   submit_script <- gsub(rm.path.prefix,'', submit_script)

   writeString( submit_script, paste0(parent.dir,'/submit_jobs.sh') )

   ## Make gather perfs script --------------------------------
   gather_script <- sprintf('#!/bin/bash
parent_dir=%s

echo "Checking if all jobs are done..."
checkAllJobsDone (){
	out_dir=$1
	done_file_path=$2

	cd $out_dir
	n_dirs=$(ls -d k*/ | wc -l)
	n_done=$(ls k*/*.job.done | wc -l)

   if [[ ! -f $done_file_path ]]; then
   	if [[ $n_dirs -eq $n_done ]]; then
   		touch $done_file_path
   	else
   		echo "Error: no. of dirs ($n_dirs) does not equal no. of done files ($n_done)"
   		exit 1
   	fi
   else
      echo "Done file exists: $done_file_path"
   fi
}
checkAllJobsDone $parent_dir/output/ $parent_dir/all.done

echo "Merging clusters.txt..."
mergeClusterTxts (){
	out_dir=$1
	out_txt=$2

	if [[ ! -f ${out_txt}.done ]]; then
		cd $out_dir

		if [[ -f $out_txt ]]; then rm $out_txt; fi

		for i in */clusters.txt.gz; do
			iter_name=$(basename $(dirname $i))
			gunzip -c $i | awk -v iter_name="$iter_name" \'{print iter_name"\t"$0}\' >> $out_txt
		done && touch ${out_txt}.done
	fi
}
mergeClusterTxts $parent_dir/output/ $parent_dir/clusters.txt.gz

echo "Merging perf.txt..."
mergePerfTxts () {
	out_dir=$1
	out_txt=$2

	if [[ ! -f ${out_txt}.done ]]; then
		cd $out_dir

		if [[ -f $out_txt ]]; then rm $out_txt; fi

		file1=$(ls */perf.txt | head -n1)
		head -n1 $file1 | gzip -c > $out_txt

		for i in */perf.txt; do
			tail -n+2 $i | gzip -c >> $out_txt
		done && touch ${out_txt}.done
	fi
}
mergePerfTxts $parent_dir/output/ $parent_dir/perfs.txt.gz
',
   sub(rm.path.prefix,'', parent.dir)
)
   writeString( gather_script, paste0(parent.dir,'/gather_perf.sh') )

}

