####################################################################################################
BIALLELE_COLS <- list(
   
   common=list(
      #sample='none',
      ensembl_gene_id='none',
      hgnc_symbol='none'
   ),
   
   mut_profile_origin=list(
      diplotype_origin='none'#,
      # a1.origin='none',
      # a2.origin='none'
   ),
   
   allele1=list(
      a1.chrom='none',
      a1.pos=0,
      a1.hgvs_c='none',
      a1.eff='none',
      a1.max_score=0,
      a1.max_score_origin='none'
   ),
   
   allele2=list(
      a2.chrom='none',
      a2.pos=0,
      a2.hgvs_c='none',
      a2.eff='none',
      a2.max_score=0,
      a2.max_score_origin='none'
   )
)

getBialleleCols <- function(groups=NULL, as=NULL){
   l <- 
      if(is.null(groups)){ unname(BIALLELE_COLS) }
      else { unname(BIALLELE_COLS[groups]) }
   
   l <- unlist(l, recursive=F)
   
   ## Return
   if(as=='names'){ names(l) } 
   else if(as=='data.frame'){ as.data.frame(l) } 
   else { l }
}

####################################################################################################
#' Merge cnv and germ/som mut profiles into a table the describes the events of allele 1 and 2
#'
#' @param mut.profile.cnv CNV mut profile table
#' @param mut.profile.mut germ or som mut profile table
#' @param mut.origin Indicate whether allele 2 originates from 'germ' or 'som'. Used to create
#' the diplotype_origin column
#' @param verbose Print progress messages?
#'
#' @return A table the describes the events of allele 1 and 2
#' @export
mkGeneDiplotypesCnvMut <- function(mut.profile.cnv, mut.profile.mut, mut.origin, verbose=T){
   # mut.profile.cnv=mut_profile$gene_cnv
   # mut.profile.mut=mut_profile$som
   # mut.origin='som'
   # mut.profile.mut=mut_profile$germ
   # mut.origin='germ'
   
   if(!(mut.origin %in% c('germ','som'))){ stop("mut.origin must be 'germ' or 'som'") }
   
   genes_cnv <- if(is.data.frame(mut.profile.cnv)){ mut.profile.cnv$ensembl_gene_id } else { NULL }
   genes_mut <- if(is.data.frame(mut.profile.mut)){ mut.profile.mut$ensembl_gene_id } else { NULL }
   
   union_genes <- unique(c(genes_cnv, genes_mut))
   n_genes <- length(union_genes)
   
   # if(verbose & n_genes!=0){
   #    pb <- txtProgressBar(min=0, max=n_genes, initial=0, style=3, width=100)
   #    counter <- 0
   # }
   
   l <- lapply(union_genes, function(i){
      # if(verbose & n_genes!=0){
      #    counter <<- counter+1
      #    setTxtProgressBar(pb,counter)
      # }
      #print(i)
      #i="ENSG00000225178"
      #i='ENSG00000012048' ## BRCA1
      #i='ENSG00000139618' ## BRCA2
      #i='ENSG00000141510' ## TP53
      #i='ENSG00000067606' ## PRKCZ
      #i='ENSG00000175279' ## CENPS
      #i='ENSG00000263001'
      
      df_cnv <- if(is.data.frame(mut.profile.cnv)){
         mut.profile.cnv[mut.profile.cnv$ensembl_gene_id==i,]
      } else {
         data.frame()
      }
      df_mut <- if(is.data.frame(mut.profile.mut)){
         mut.profile.mut[mut.profile.mut$ensembl_gene_id==i,]
      } else {
         data.frame()
      }
      
      if(nrow(df_cnv)==0 & nrow(df_mut)==0){
         #return( getBialleleCols(as='data.frame') )
         return( NULL )
      }
      
      #========= Full gene loss / truncation cases =========#
      if( nrow(df_cnv)!=0 && df_cnv$loss_type=='deep_deletion' ){
         ## df_cnv should only ever have one row
         out <- getBialleleCols(as='data.frame')
         
         ## Scores
         out$a1.max_score <- df_cnv$a1.score
         out$a2.max_score <- df_cnv$a2.score
         
         ## Variant effect
         out$a1.eff <- out$a2.eff <- df_cnv$loss_type
         
         ## Chrom / pos
         out$a1.chrom <- out$a2.chrom <- df_cnv$chrom
         out$a1.pos <- df_cnv$start
         out$a2.pos <- df_cnv$end
         
         out <- within(out,{
            a1.max_score_origin <- a2.max_score_origin <- 'cnv'
            #a1.origin <- a2.origin <- 'cnv'
            diplotype_origin <- 'cnv_cnv'
         })
         
         ## Common
         common_cols <- getBialleleCols('common',as='names')
         out[common_cols] <- df_cnv[common_cols]
         
         return(out)
      }
      
      #========= Main =========#
      common_cols <- getBialleleCols('common',as='names')
      
      #--------- CNV (left hand side) ---------#
      out_cnv <- getBialleleCols(c('common','mut_profile_origin','coords','allele1'),as='data.frame')
      
      if(nrow(df_cnv)!=0){
         out_cnv$a1.chrom <- df_cnv$chrom
         out_cnv$a1.max_score <- df_cnv$a1.score
         out_cnv$a1.eff <- df_cnv$loss_type
         out_cnv$a1.max_score_origin <- 'cnv'
      }
      
      ## fill origin cols
      out_cnv$diplotype_origin <- paste0('cnv_',mut.origin)
      # out_cnv <- within(out_cnv, {
      #    diplotype_origin <- paste0('cnv_',mut.origin)
      #    a1.origin <- 'cnv'
      #    a2.origin <- mut.origin
      # })
      
      #--------- germ or som (right hand side) ---------#
      out_mut <- getBialleleCols(c('common','coords','allele2'),as='data.frame')
      
      if(nrow(df_mut)!=0){
         ## Initialize empty df
         out_mut <- getBialleleCols(c('common','coords','allele2'),as='data.frame')
         out_mut <- out_mut[rep(1,nrow(df_mut)),]
         out_mut[common_cols] <- df_mut[common_cols]
         
         ## Fill metadata
         meta_cols <- getBialleleCols(c('common','coords'),as='names')
         out_mut[meta_cols] <- df_mut[meta_cols]
         
         ## Fill allele data
         allele_cols_out <- getBialleleCols('allele2',as='names')
         
         allele_cols_source <- allele_cols_out
         allele_cols_source[allele_cols_source=='a2.eff'] <- 'snpeff_eff'
         allele_cols_source <- gsub('a2[.]','',allele_cols_source)
         
         out_mut[allele_cols_out] <- df_mut[allele_cols_source]
      }
      
      ## Add gene names
      out_cnv[common_cols] <- if(nrow(df_cnv)!=0){
         df_cnv[common_cols]
      } else {
         df_mut[1,common_cols]
      }
      
      #--------- Join ---------#
      #out <- merge(out_cnv, out_mut, all=T)
      out <- cbind(out_cnv, out_mut[!(colnames(out_mut) %in% common_cols)])
      
      ## column check
      #ncol(out)
      #length(getBialleleCols(as='names'))
      
      return(out)
   })
   
   return( do.call(rbind,l) )
}

# df_cnv_mut <- mergeCnvMut(mut_profile$cnv,mut_profile$germ,'germ')
# View(df_cnv_mut)
# table(df_cnv_mut$a1)

####################################################################################################
#' Merge germ/som mut profiles into a table the describes the events of allele 1 and 2
#'
#' @description Per gene, only the most pathogenic germ/som variant pairs are merged. This
#' prevents the table from growing exponentially.
#'
#' @param mut.profile.mut1 germ/som mut profile table
#' @param mut.profile.mut2 germ/som mut profile table
#' @param diplotype.origin Indicate the mut/mut combination: 'germ_som', 'som_som', 'germ_germ'
#' @param verbose Print progress messages?
#'
#' @return A table the describes the events of allele 1 and 2
#' @export
mkGeneDiplotypesMutMut <- function(
   mut.profile.mut1, mut.profile.mut2, 
   diplotype.origin, min.biall.score=c(3,3),
   get.max.effect.per.gene=T,
   verbose=T
){
   #mut.profile.mut1=mut_profile$som
   #mut.profile.mut1=mut_profile$germ
   #mut.profile.mut2=mut_profile$som
   
   if(nrow(mut.profile.mut1)==0 | nrow(mut.profile.mut2)==0){
      return(NULL)
   }
   
   #========= Col names =========#
   common_cols <- getBialleleCols('common',as='names')
   
   sel_cols <- list(
      common=getBialleleCols('common',as='names'),
      a1_out=getBialleleCols('allele1',as='names'),
      a2_out=getBialleleCols('allele2',as='names')
   )
   
   sel_cols$a1_source <- (function(){
      v <- sel_cols$a1_out
      v[v=='a1.eff'] <- 'snpeff_eff'
      v <- gsub('a1[.]','',v)
      return(v)
   })()
   
   sel_cols$a2_source <- (function(){
      v <- sel_cols$a2_out
      v[v=='a2.eff'] <- 'snpeff_eff'
      v <- gsub('a2[.]','',v)
      return(v)
   })()
   
   #========= Merge =========#
   formatMutProfile <- function(mut.profile.mut, allele){
      df <- mut.profile.mut[, unlist(sel_cols[ c('common',paste0(allele,'_source') )], use.names=F) ]
      
      if(allele=='a1'){
         df <- df[df$max_score>=min.biall.score[1],]
      } else if(allele=='a2') {
         df <- df[df$max_score>=min.biall.score[2],]
      }
      
      colnames(df) <- unlist(sel_cols[c('common',paste0(allele,'_out'))], use.names=F)
      return(df)
   }
   
   df_mut1 <- if(is.data.frame(mut.profile.mut1)){
      formatMutProfile(mut.profile.mut1,'a1')
   } else {
      getBialleleCols(c('common','allele1'),as='data.frame')
   }
   
   df_mut2 <- if(is.data.frame(mut.profile.mut2)){
      formatMutProfile(mut.profile.mut2,'a2')
   } else {
      getBialleleCols(c('common','allele2'),as='data.frame')
   }
   
   out <- merge(df_mut1, df_mut2, by=sel_cols$common, all=T)
   
   if(nrow(out)==0){
      return(NULL)
   }
   
   #========= Fill NA's =========#
   na_fills <- getBialleleCols(c('common','allele1','allele2'), as='data.frame')
   
   for(i in colnames(na_fills)){
      out[is.na(out[,i]),i] <- na_fills[,i]
   }
   
   #========= Remove dup rows =========#
   if(identical(mut.profile.mut1,mut.profile.mut2)){
      mkVariantStrings <- function(df){
         lapply(c('a1','a2'),function(i){
            #i='a1'
            cols <- c('ensembl_gene_id', paste0(i,'.hgvs_c'))
            apply(df[cols], 1, function(j){
               paste0(j, collapse=';')
            })
         })
      }
      
      ## Remove duplicate rows where the two variants are exactly the same
      variant_strings <- mkVariantStrings(out)
      out <- out[variant_strings[[1]] != variant_strings[[2]],]
      
      if(nrow(out)==0){
         return(NULL)
      }
      
      ## Remove duplicate rows where a1 and a2 are mirrored
      variant_strings <- mkVariantStrings(out)
      biall_strings1 <- paste0(variant_strings[[1]],';',variant_strings[[2]])
      biall_strings2 <- paste0(variant_strings[[2]],';',variant_strings[[1]])
      
      out$biall_strings1 <- biall_strings1
      
      i<-1
      while(i<length(biall_strings1)){
         dup_index <- which(biall_strings2 %in% biall_strings1[i])
         biall_strings1 <- biall_strings1[-dup_index]
         biall_strings2 <- biall_strings2[-dup_index]
         i <- i+1
      }
      
      out <- out[out$biall_strings1 %in% biall_strings1,]
      out$biall_strings1 <- NULL
   }
   
   #========= Export =========#
   ## Select max pathogenicity only
   if(get.max.effect.per.gene){
      out_split <- split(out, out$hgnc_symbol)
      out <- do.call(rbind,lapply(out_split,function(i){
         max_index <- which.max(rowSums(i[,c('a1.max_score','a2.max_score')]))
         i[max_index,]
      }))
   }
   
   ## Add diplotype origin column
   out <- insColAfter(out, diplotype.origin, after='hgnc_symbol',colname='diplotype_origin')
   
   return(out)
}




