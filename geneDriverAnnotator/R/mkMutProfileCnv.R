#' Determine gene copy number gains and losses
#'
#' @description Annotations deep deletions, truncations, LOH, and amplifications
#' 
#' @param cnv.file Path to the cnv bed file.
#' @param sel.cols A named vector indicating the columns to keep. Destination colnames should be:
#' chrom, start, end, total_cn, major_cn, minor_cn
#' @param exons.bed.file Path to the bed file containing the columns chrom, exon_start, exon_end, 
#' hgnc_symbol, ensembl_gene_id
#' @param genes.bed.file Alternative to exons.bed.file.
#' @param min.gain.ratio.genome The min genome amp ratio (min_copy_number / genome_ploidy) to be 
#' considered an amplification.
#' @param min.gain.ratio.max.arm.diff The min difference between gain_ratio_max and gain_ratio_arm
#' to be considered a focal amplification
#' @param deep.del.max.min.copy.number The max min_copy_number for a gene to be considered to 
#' be completely lost (deep deletion)
#' @param trunc.max.min.copy.number The max min_copy_number to for a gene to be considered to have
#' a truncation
#' @param loh.max.min.minor.allele.ploidy The maximum min_minor_allele_ploidy for a gene to be
#' considered to have loss of heterozygosity
#'
#' @return The original input dataframe with annotation columns
#' @export
#'
mkMutProfileGeneCnv <- function(
   cnv.file=NULL, 
   #gene.cnv.file=NULL,
   sel.cols=c(chrom='chromosome',start='start',end='end',total_cn='copyNumber',major_cn='majorAllelePloidy',minor_cn='minorAllelePloidy'),
   exons.bed.file=EXONS_BED_FILE, genes.bed.file=NULL, 
   
   ## Cutoffs
   min.gain.ratio.genome=1.8,
   min.gain.ratio.max.arm.diff=0.8,
   
   deep.del.max.min.copy.number=0.3, 
   loh.max.min.minor.allele.ploidy=0.2,
   
   ## Legacy args
   #genes.enst2ensg.file=GENES_ENST2ENSG,
   
   ## Misc
   verbose=T
){
   ## ZNF703 focal amp:
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/170907_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv'
   ## MYC arm gain:
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/190115_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv'
   ## ZNF217 chrom gain:
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/180616_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv'
   ## No CNA:
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data//somatics/170613_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv'
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/161123_HMFregXXXXXXXX/XXXXXXXX.purple.cnv.somatic.tsv'
   
   #cnv.file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/PCAWG_2020/vcf/somatic/cnv/cna_annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt'
   #cnv.file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/PCAWG_2020/vcf/somatic/cnv/cna_annotated/fc8130df-6860-7677-e040-11ac0d485ddc.consensus.20170119.somatic.cna.annotated.txt'
   #sel.cols=c(chrom='chromosome',start='start',end='end',total_cn='total_cn',major_cn='major_cn',minor_cn='minor_cn')
   #cnv.file<-input.file.paths['cnv']
   
   #cnv.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104-update3//somatics/181019_HMFregXXXXXXXX/purple/XXXXXXXX.purple.cnv.somatic.tsv'
   #sel.cols=c(chrom='chromosome',start='start',end='end',total_cn='copyNumber',major_cn='majorAlleleCopyNumber',minor_cn='minorAlleleCopyNumber')
   #exons.bed.file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/Chromatin_modifiers/gene_list/chromatin_modifier_exons.bed.gz'
   
   ## ------------------------
   if(verbose){ message('Reading bed file...') }
   if(!is.null(exons.bed.file)){
      bed <- read.delim(exons.bed.file, stringsAsFactors=F)
      colnames(bed)[2:3] <- c('start','end')
      gene_info <- bed[,c('chrom','gene_start','gene_end','hgnc_symbol','ensembl_gene_id')]
      colnames(gene_info)[2:3] <- c('start','end')
   } else {
      bed <- read.delim(genes.bed.file, stringsAsFactors=F)
      gene_info <- bed[,c('chrom','start','end','hgnc_symbol','ensembl_gene_id')]
   }
   gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id),]
   gene_info$chrom_arm <- getChromArm(gene_info[,c('chrom','start','end')], seq.levels.style='NCBI')
   
   ## Subsetting for genes ------------------------
   #if(is.null(gene.cnv.file)){
   if(verbose){ message('Reading cnv file...') }
   cnv <- read.delim(cnv.file, stringsAsFactors=F, check.names=F)
   colnames(cnv) <- sub('#','', colnames(cnv))

   cnv <- selectRequiredCols(
      df=cnv,
      required.cols=c('chrom','start','end','total_cn','major_cn','minor_cn'),
      sel.cols=sel.cols
   )
   
   GenomeInfoDb::seqlevelsStyle(cnv$chrom)<- 'NCBI'
   
   if(verbose){ message('Finding CNVs that overlap with bed file exons...') }
   overlaps <- isOverlappingChromPos(
      df1=bed[,c('chrom','start','end')],
      df2=cnv[,c('chrom','start','end')],
   )
   
   if(verbose){ message('Calculating CN info per gene...') }
   
   if(any(sapply(overlaps,length)==0)){
      warning('Some genes have no copy number info')
   }
   
   if(!is.null(exons.bed.file)){
      exon_cnv <- lapply(seq_along(overlaps),function(i){
         #i=1
         if(length(overlaps[[i]])==0){
            return(NULL)
         }
         df <- cnv[overlaps[[i]],]
         df$hgnc_symbol <- bed[i,'hgnc_symbol']
         df$ensembl_gene_id <- bed[i,'ensembl_gene_id']
         return(df)
      })
      exon_cnv <- do.call(rbind, exon_cnv)
      
      exon_cnv_split <- split(exon_cnv, factor(exon_cnv$ensembl_gene_id, unique(exon_cnv$ensembl_gene_id)))
      suppressWarnings({
         gene_cnv <- do.call(rbind, lapply(exon_cnv_split, function(i){
            #i=exon_cnv_split[['ENSG00000168172']]
            data.frame(
               ensembl_gene_id=i[1,'ensembl_gene_id'],
               min_copy_number=min(i$total_cn, na.rm=T),
               max_copy_number=max(i$total_cn, na.rm=T),
               min_minor_allele_ploidy=min(i$minor_cn, na.rm=T)
            )
         }))
      })
      
   } else {
      gene_cnv_pre <- lapply(seq_along(overlaps),function(i){
         #i=1
         if(length(overlaps[[i]])==0){
            return(NULL)
         }
         df <- cnv[overlaps[[i]],]
         df$hgnc_symbol <- bed[i,'hgnc_symbol']
         df$ensembl_gene_id <- bed[i,'ensembl_gene_id']
         return(df)
      })
      gene_cnv_pre <- do.call(rbind, gene_cnv_pre)
      gene_cnv_split <- split(gene_cnv_pre, factor(gene_cnv_pre$ensembl_gene_id), unique(gene_cnv_pre$ensembl_gene_id))
      gene_cnv <- do.call(rbind, lapply(gene_cnv_split, function(i){
         data.frame(
            ensembl_gene_id=i[1,'ensembl_gene_id'],
            min_copy_number=min(i$total_cn, na.rm=T),
            max_copy_number=max(i$total_cn, na.rm=T),
            min_minor_allele_ploidy=min(i$minor_cn, na.rm=T)
         )
      }))
   }
   rownames(gene_cnv) <- NULL
   
   gene_cnv$min_copy_number[!is.finite(gene_cnv$min_copy_number)] <- NA
   gene_cnv$max_copy_number[!is.finite(gene_cnv$max_copy_number)] <- NA
   gene_cnv$min_minor_allele_ploidy[!is.finite(gene_cnv$min_minor_allele_ploidy)] <- NA
   
   gene_cnv <- cbind(
      gene_info,
      gene_cnv[match(gene_info$ensembl_gene_id, gene_cnv$ensembl_gene_id),-1]
   )
   # } else {
   #    gene_cnv_pre <- read.delim(gene.cnv.file, check.names=F)
   #    
   #    ## Format columns
   #    sel_cols=c(
   #       chrom='Chromosome',start='Start',end='End',
   #       min_copy_number='MinCopyNumber',max_copy_number='MaxCopyNumber',
   #       min_minor_allele_ploidy='MinMinorAllelePloidy',
   #       ensembl_transcript_id='TranscriptId'
   #    )
   #    gene_cnv_pre <- gene_cnv_pre[sel_cols]
   #    colnames(gene_cnv_pre) <- names(sel_cols)
   #    
   #    ## Get ENSG ids
   #    genes_enst2ensg <- read.delim(GENES_ENST2ENSG)
   #    ins_cols <- genes_enst2ensg[
   #       match(gene_cnv_pre$ensembl_transcript_id, genes_enst2ensg$ensembl_transcript_id),
   #       c('ensembl_gene_id','hgnc_symbol')
   #    ]
   #    gene_cnv <- cbind(gene_cnv_pre,ins_cols)
   #    gene_cnv <- gene_cnv[gene_cnv$ensembl_gene_id %in% bed$ensembl_gene_id,]
   #    gene_cnv <- gene_cnv[order(gene_cnv$hgnc_symbol),]
   #    
   #    ## Load CNV 
   #    cnv <- read.delim(cnv.file, stringsAsFactors=F, check.names=F)
   #    colnames(cnv) <- sub('#','', colnames(cnv))
   #    cnv <- cnv[,sel.cols]
   #    colnames(cnv) <- names(sel.cols)
   #    GenomeInfoDb::seqlevelsStyle(cnv$chrom)<- 'NCBI'
   # }
   
   ## ----------------------------------------------------------------
   if(verbose){ message('Calculating arm ploidies...') }
   sel_cols <- structure(colnames(cnv), names=colnames(cnv)) 
   arm_ploidies <- data.frame(
      total_cn=calcChromArmPloidies(cnv=cnv, mode='total_cn', sel.cols=sel_cols),
      minor_cn=calcChromArmPloidies(cnv=cnv, mode='minor_cn', sel.cols=sel_cols)
   )
   genome_ploidy <- unlist(arm_ploidies['genome',])
   arm_ploidies <- arm_ploidies[rownames(arm_ploidies) != 'genome',]
   
   if(verbose){ message('Identifying whole chromosome CNAs...') }
   whole_chrom_cna <- (function(){
      chrom_names <- gsub('p|q','',rownames(arm_ploidies))
      chrom_names <- factor(chrom_names,unique(chrom_names))
      chrom_gain <- unlist(lapply(split(arm_ploidies$total_cn, chrom_names), function(i){
         if(length(i)==1 & i[1] > genome_ploidy['total_cn']){ return(TRUE) }
         if( i[1]==i[2] & i[1] > genome_ploidy['total_cn'] ){ return(TRUE) }
         return(FALSE)
      }))
      chrom_loh <- unlist(lapply(split(arm_ploidies$minor_cn, chrom_names), function(i){
         if(length(i)==1 & i[1]==0){ return(TRUE) }
         if( i[1]==0 & i[2]==0 ){ return(TRUE) }
         return(FALSE)
      }))
      
      data.frame(chrom_gain, chrom_loh)
   })()
   
   ## ----------------------------------------------------------------
   if(verbose){ message('> Gains...') }
   gains <- gene_cnv
   gains$genome_ploidy <- genome_ploidy[['total_cn']]
   gains$arm_ploidy <- arm_ploidies[gene_cnv$chrom_arm,'total_cn']
   
   if(verbose){ message('Calculating amplification levels...') }
   gains$gain_ratio_genome <- gains$min_copy_number / gains$genome_ploidy
   gains$gain_ratio_arm <- gains$arm_ploidy / gains$genome_ploidy
   gains$gain_ratio_focal <- gains$min_copy_number / gains$arm_ploidy
   
   gains$gain_ratio_max <- pmax(gains$gain_ratio_genome, gains$gain_ratio_arm, gains$gain_ratio_focal)
   
   
   if(verbose){ message('Determining amp type...') }
   gains$gain_type <- ifelse(
      gains$gain_ratio_genome >= min.gain.ratio.genome,
      'amp','none'
   )
   gains$gain_type[ is.na(gains$gain_type) ] <- 'none'
   
   gains$gain_type[
      gains$gain_type=='amp' & 
      (gains$gain_ratio_max - gains$gain_ratio_arm) < min.gain.ratio.max.arm.diff
   ] <- 'arm'
   
   gains$gain_type[
      gains$gain_type=='arm' & 
      whole_chrom_cna[match(gains$chrom, rownames(whole_chrom_cna)),'chrom_gain']
   ] <- 'chrom'
   
   gains$gain_type[
      gains$gain_type=='amp'
   ] <- 'focal'
   
   # if(F){
   #    linx_genes <- read.delim(
   #       '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/hmfSvFeatureExtractor/scripts/gather_linx_driver_tsvs/linx_driver_genes.txt',
   #       header=F, stringsAsFactors=F
   #    )[,1]
   #    
   #    gains[order(gains$gain_ratio_max, decreasing=T),]
   #    gains[gains$hgnc_symbol %in% linx_genes & gains$gain_type!='none',]
   # }

   ## ----------------------------------------------------------------
   if(verbose){ message('> Losses...') }
   losses <- gene_cnv
   
   if(verbose){ message('Determining loss type...') }
   losses$loss_type <- 'none'
   losses <- within(losses,{
      loss_type[ min_minor_allele_ploidy <= loh.max.min.minor.allele.ploidy ] <- 'loh'
      loss_type[ min_copy_number <= deep.del.max.min.copy.number ] <- 'deep_deletion'
   })
   
   if(verbose){ message('Calculating biallelic loss scores...') }
   loss_scores <- as.data.frame(matrix(
      0, nrow=nrow(losses), ncol=2, 
      dimnames=list(NULL, c('a1.score','a2.score'))
   ))
   loss_scores[losses$loss_type %in% c('loh','deep_deletion'),'a1.score'] <- 5
   loss_scores[losses$loss_type=='deep_deletion','a2.score'] <- 5
   
   if(verbose){ message('Determining detailed loss types...') }
   losses$arm_minor_cn <- arm_ploidies[match(losses$chrom_arm, rownames(arm_ploidies)),'minor_cn']
   
   losses$loss_type[
      losses$loss_type=='loh' & 
      losses$arm_minor_cn==0
   ] <- 'loh_arm'
   
   losses$loss_type[
      losses$loss_type %in% c('loh','loh_arm') & 
      whole_chrom_cna[match(losses$chrom, rownames(whole_chrom_cna)),'chrom_loh']
   ] <- 'loh_chrom'
   
   losses <- cbind(losses, loss_scores)
   
   out <- cbind(
      gene_cnv,
      gains[!(colnames(gains) %in% colnames(gene_cnv))],
      losses[!(colnames(losses) %in% colnames(gene_cnv))]
   )
   
   ## Round floats to save space
   out <- as.data.frame(lapply(out, function(i){
      if(is.numeric(i)){ round(i,3) }
      else { i }
   }))

   return(out)
}





