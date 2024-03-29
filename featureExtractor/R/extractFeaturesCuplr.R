#' Extract input features for CUPLR
#'
#' @param in.dir A directory containing the input files. The required files are searched for with
#' the following regular expressions:
#' \itemize{
#' \item  purple.smnv:        '.purple.somatic.vcf.gz$'
#' \item  purple.cnv:         '.purple.cnv.somatic.tsv$'
#' \item  purple.purity:      '.purple.purity.tsv$'
#' \item  linx.drivers:       '.linx.driver.catalog.tsv$'
#' \item  linx.fusions:       '.linx.fusion.tsv$'
#' \item  linx.viral.inserts: '.linx.viral_inserts.tsv$'
#' \item  linx.vis.sv.data:   '.linx.vis_sv_data.tsv$'
#' }
#' @param input.paths A named list of the input paths. See `in.dir` for the required names 
#' @param out.dir (Optional) A directory to write intermediate and final output files. If
#' unspecified, a dataframe of the features will be returned
#' @param clonal.variants.only If TRUE, only clonal variants will be used for the mutational 
#' signature and RMD signature features
#' @param verbose Show progress messages? 0: No messages. 1: Messages from `extractFeaturesCuplr()`. 
#' 2: Messages from internal functions called by `extractFeaturesCuplr()`
#'
#' @return See `out.dir`
#' @export
#'
extractFeaturesCuplr <- function(
   in.dir=NULL, input.paths=NULL, out.dir=NULL, clonal.variants.only=F, verbose=F
){
   
   ## Debugging --------------------------------
   if(F){
      devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/featureExtractor')
      in.dir='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/doc/data/DO35886/'
      out.dir=paste0(in.dir,'/output/')
      dir.create(out.dir, showWarnings=F)
      verbose=T
   }
   
   #print(input.paths)
   
   ## Init ================================
   ## Inputs --------------------------------
   input_path_patterns <- c(
      purple.smnv=        '.purple.somatic.vcf.gz$',
      purple.cnv=         '.purple.cnv.somatic.tsv$',
      purple.purity=      '.purple.purity.tsv$',
      linx.drivers=       '.linx.driver.catalog.tsv$',
      linx.fusions=       '.linx.fusion.tsv$',
      linx.viral.inserts= '.linx.viral_inserts.tsv$',
      linx.vis.sv.data=   '.linx.vis_sv_data.tsv$'
   )
   
   if(!is.null(in.dir)){
      if(!dir.exists(in.dir)){ stop("`in.dir` does not exist: ", in.dir) }
      input.paths <- sapply(input_path_patterns, function(i){
         list.files(path=in.dir, pattern=i, full.names=T)
      })
   }
   
   if(!is.null(input.paths)){
      
      input.paths[input.paths=="NA"] <- NA
      
      if(!all(names(input.paths) %in% names(input_path_patterns))){
         stop(
            "`input.paths` must have the names: ", 
            paste(names(input_path_patterns), collapse=', ')
         )
      }
      
      valid_input_paths <- sapply(input.paths, function(i){ is.na(i) || length(i)!=0 || file.exists(i) })
      if(!all(valid_input_paths)){
         stop(
            "Some `input.paths` do not exist:\n", 
            paste(input.paths[!valid_input_paths], collapse='\n')
         )
      }
      
   } else {
      stop('Input must be specified to `input.paths` or `in.dir`')
   }
   
   ## Helper functions --------------------------------
   ## Code for saving progress
   if(!is.null(out.dir)){
      if(!dir.exists(out.dir)){ stop('`out.dir` does not exist') }
      setwd(out.dir)
      dir.create('raw', showWarnings=F)
      dir.create('features', showWarnings=F)
   }
   saveAndReadVector <- function(v, path){
      #v=contexts
      #path='/Users/lnguyen/Desktop/test.txt'
      
      if(is.null(out.dir)){ return(v) }
      
      if(!file.exists(path)){
         if(verbose>=2){ message('Writing output to: ', path) }
         write.table(
            x=matrix(v, ncol=1, dimnames=list(names(v), NULL)),
            file=path, sep='\t', col.names=F, quote=F
         )
         return(v)
      }
      
      if(verbose>=2){ message('Reading output from: ', path) }
      df <- read.delim(path, header=F, check.names=F)
      structure(df[,2], names=df[,1])
   }
   
   ##
   normalizeVector <- function(v){
      v <- v / sum(v)
      v[is.na(v)] <- 0
      return(v)
   }
   
   ## Extract each feature type ================================
   features <- list()
   
   ## Load the SMNV vcf once (used for multiple feature types)
   vcf_smnv <- mutSigExtractor::variantsFromVcf(
      vcf.file=input.paths[['purple.smnv']], 
      vcf.fields=c(1,2,4,5,7,8),
      vcf.filter='PASS', keep.chroms=c(1:22,'X'),
      verbose=verbose>=2
   )
   
   ## --------------------------------
   if(verbose){ message('> SNV/indel/DBS signatures') }
   contexts <- saveAndReadVector(
      extractContextsSmnvIndel(df=vcf_smnv, as.matrix=F, clonal.variants.only=clonal.variants.only, verbose=verbose>=2),
      'raw/smnv_contexts.txt'
   )
   
   contexts_split <- splitFeaturesByGroup.default(contexts, rm.tags=T)
  
   sigs <- list(
      snv=mutSigExtractor::fitToSignatures(
         mut.context.counts=contexts_split$snv, 
         signature.profiles=featureExtractor::SBS_SIGNATURE_PROFILES_V3
      ),
      
      dbs=mutSigExtractor::fitToSignatures(
         mut.context.counts=contexts_split$dbs,
         signature.profiles=featureExtractor::DBS_SIGNATURE_PROFILES
      ),
      
      indel=mutSigExtractor::fitToSignatures(
         mut.context.counts=contexts_split$indel,
         signature.profiles=featureExtractor::INDEL_SIGNATURE_PROFILES
      )
   )

   sigs$snv <- combineFeatures.numeric(sigs$snv, regex='^SBS7', target.name='SBS7')
   sigs$snv <- combineFeatures.numeric(sigs$snv, regex='^SBS10', target.name='SBS10')
   sigs$snv <- combineFeatures.numeric(sigs$snv, regex='^SBS17', target.name='SBS17')
   
   sigs <- lapply(sigs, normalizeVector)
   
   features$sigs <- saveAndReadVector(
      do.call(c, unname(sigs)),
      'features/sigs.txt'
   )
   rm(sigs)
   
   features$mut_load <- saveAndReadVector(
      sapply(contexts_split, sum),
      'features/mut_load.txt'
   )
   
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> Regional mutational density') }
   rmd_bin_counts <- saveAndReadVector(
      extractRmd(df=vcf_smnv, clonal.variants.only=clonal.variants.only, as.matrix=F, verbose=verbose>=2),
      'raw/rmd_bin_counts.txt'
   )
   
   #features$rmd <- rmd_bin_counts
   features$rmd <- saveAndReadVector(
      normalizeVector(rmd_bin_counts),
      'features/rmd.txt'
   )
   
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> SV contexts') }
   features$sv <- saveAndReadVector(
      extractContextsSvLinx(input.paths[['linx.vis.sv.data']]),
      'features/sv.txt'
   )
   
   ## --------------------------------
   if(verbose){ message('> Chrom arm ploidies') }
   arm_ploidy <- saveAndReadVector(
      calcChromArmPloidies(
         cnv.file=input.paths[['purple.cnv']],
         sel.cols=c(chrom='chromosome',start='start',end='end',total_cn='copyNumber',major_cn='majorAlleleCopyNumber',minor_cn='minorAlleleCopyNumber'),
         verbose=verbose>=2
      ),
      'raw/arm_ploidy.txt'
   )
      
   chrom_arm <- calcChromArmCnChange(arm_ploidy, direction='both')
   
   features$chrom_arm <- saveAndReadVector(
      unlist(chrom_arm),
      'features/chrom_arm.txt'
   )
   
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> Gene drivers') }
   features$gene <- saveAndReadVector(
      getGeneDriverEvents(input.paths[['linx.drivers']]),
      'features/gene.txt'
   )
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> Gene fusions') }
   features$fusion <- saveAndReadVector(
      getGeneFusions(input.paths[['linx.fusions']]),
      'features/fusion.txt'
   )
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> Viral insertions') }
   features$viral_ins <- saveAndReadVector(
      getViralInsertions(input.paths[['linx.viral.inserts']]),
      'features/viral_ins.txt'
   )
   if(verbose>=2){ message() }
   
   ## --------------------------------
   if(verbose){ message('> PURPLE purity') }
   purple <- saveAndReadVector(
      getPurplePurityData(input.paths[['purple.purity']]),
      'features/purple.txt'
   )
   
   features$genome <- purple[c('genome_ploidy','diploid_proportion','has_wgd')]
   names(features$genome)[names(features$genome)=='genome_ploidy'] <- 'ploidy'
   
   features$gender <- purple['gender']
   
   if(verbose>=2){ message() }
   
   ## Output ================================
   features <- lapply(features, function(i){
      #i=features[[1]]
      as.data.frame(
         matrix(i, nrow=1, dimnames=list(NULL,names(i)) ),
         check.names=F
      )
   })
   features <- as.data.frame(
      unlist(features, recursive=F), 
      check.names=F
   )
   
   ## Fix bool features that were coerced to numeric due to being in a vector
   features$gender.gender <- as.logical(features$gender.gender)
   features$genome.has_wgd <- as.logical(features$genome.has_wgd)
   
   if(!is.null(out.dir)){
      if(verbose){ message('> Features saved to: ',out.dir,'/all_features.txt.gz') }
      write.table(features, gzfile('all_features.txt.gz'), sep='\t', quote=F, row.names=F)
      return(invisible(NULL))
   }
   
   return(features)
}
