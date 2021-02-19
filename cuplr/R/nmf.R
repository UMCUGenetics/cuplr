#' Non-negative matrix factorization
#'
#' @description A wrapper function around NNLM::nnmf() to extract signatures from a numeric matrix
#' (where columns are features and rows are samples).
#'
#' The function determines the optimum NMF rank by first setting a portion of values
#' (`impute.prop`) in the input matrix to NA, which triggers value imputation from NNLM::nnmf(). A
#' mean squared error (MSE) between the original and imputed values is determined. The above
#' procedure is repeated several times (`repeats`), after which median MSE values are calculated.
#' The optimum rank is in theory the rank at which MSE is minimum, however in practice, this is
#' the rank before which MSE increases by a certain threshold (`max.rel.log.mse.increase`).
#'
#' @param A A numeric matrix (rows=samples, columns=features)
#' @param k.range Integer vector. A range of ranks in which to find the optimum rank
#' @param repeats See description
#' @param max.samples Subsample number of rows to this value to reduce computation time
#' @param impute.prop See description
#' @param max.rel.log.mse.increase See description
#' @param seed Random seed
#' @param verbose Show progress? Can be 0,1,2
#'
#' @return A list containing the factorized matrices and MSE values from rank search
#' @export
#'
nmf <- function(
   A,
   k.range=1:10, repeats=10, max.samples=100, impute.prop=0.1, max.rel.log.mse.increase=0.002,
   run.with.optimum.k=T, seed=1, verbose=0
){
   if(F){
      #m_rmd <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.10e_rmd1mb_DR104update/features/features.txt.gz')
      A=subset(m_rmd,response=='Breast')[,-1]
      k.range=1:10
      repeats=10
      max.samples=100
      impute.prop=0.1
      max.rel.log.mse.increase=0.002
      seed=1
      verbose=T
   }

   ## Init --------------------------------
   require(NNLM)

   set.seed(seed)

   A <- as.matrix(A)
   if(!is.numeric(A)){ stop('`A` must be a numeric matrix or dataframe') }

   ## Calc NMF performance; Do final fit --------------------------------
   main <- function(k, repeat.num, return.perf=F, max.samples){

      i_seed <- as.integer(paste0(seed,repeat.num))
      set.seed(i_seed)

      ## Randomly replace some values with NA
      ## nnmf() performs imputation
      A2 <- A

      ## Reduce number of samples to reduce run time
      if(!is.null(max.samples) & nrow(A2)>max.samples){
         A2 <- A2[sample(nrow(A2), max.samples),]
      }

      if(return.perf){
         ind <- sample(length(A2), impute.prop*length(A2))
         A2[ind] <- NA
      }

      ## Check k
      min.k <- min(dim(A2))
      A.isNA <- is.na(A2)
      A.anyNA <- any(A.isNA) # anyNA is depreciated in new version of R
      #print('1')
      if (A.anyNA) {
         min.k <- min(min.k, ncol(A2) - rowSums(A.isNA), nrow(A2) - colSums(A.isNA))
      }
      #print('2')
      #print(k)
      #print(min.k)
      if(k > min.k){
         warning(paste0('k is larger than min.k (',min.k,'). Returning NA'))
         return(
            data.frame(
               mse_imputed=NA,
               repeat_num=repeat.num
            )
         )
      }
      #print('3')

      ## Do NMF
      nnmf_out <- NNLM::nnmf(
         A=A2, k=k, loss='mkl', max.iter=2000L,
         verbose=if(verbose==2){ 2 } else { 0 }
      )

      if(return.perf){
         mse_imputed <- mean(((nnmf_out$W %*% nnmf_out$H)[ind] - A[ind])^2 )

         return(data.frame(
            mse_imputed=mse_imputed,
            repeat_num=repeat.num
         ))
      }

      return(nnmf_out)
   }

   ## Determine optimal rank --------------------------------
   if(verbose>0){
      if(length(k.range)>1){
         message('Calculating MSE and determining optimal rank...')
      } else {
         message('Calculating MSE...')
      }
   }
   param_grid <- data.frame(
      k=rep(k.range, each=repeats),
      repeat_num=rep(1:repeats, length(k.range)),
      index=1:(repeats*length(k.range))
   )

   ## At the below indexes, calculate median mse
   checkpoints <- param_grid$index[ param_grid$index %% repeats == 0 ]

   perf <- data.frame(k=numeric(), repeat_num=numeric(), mse_imputed=numeric())
   perf_summ <- data.frame(k=numeric(), log10_mse_imputed=numeric())

   first_k <- min(param_grid$k)
   for(i in 1:nrow(param_grid)){
      k <- param_grid[i,'k']
      repeat_num <- param_grid[i,'repeat_num']

      if(verbose>0){ message('> Rank ',k,', repeat ',repeat_num) }

      i_perf <- main(k=k, repeat.num=repeat_num, return.perf=T, max.samples=max.samples)
      perf <- rbind(perf, cbind(k=k,i_perf))

      if(i %in% checkpoints){
         k_perf <- perf[perf$k==k,]
         perf_summ <- rbind(
            perf_summ,
            data.frame(
               k=k,
               log10_mse_imputed=log10(median(k_perf$mse_imputed, na.rm=T))
            )
         )

         if(k==first_k){ next }

         perf_before <- perf_summ[(nrow(perf_summ)-1),'log10_mse_imputed']
         perf_after <- perf_summ[nrow(perf_summ),'log10_mse_imputed']

         mse_increase <- (perf_before - perf_after)/perf_before

         if(!is.null(max.rel.log.mse.increase)){
            if(mse_increase > max.rel.log.mse.increase){ break }
         }
      }
   }

   ## Run NMF with optimum rank --------------------------------
   ## Optimum rank: Rank before mse increases above threshold (i.e. max.rel.log.mse.increase)

   if(length(k.range)==1){
      optimum_k <- k.range
   } else {
      if(tail(k.range,1)==nrow(perf_summ)){
         optimum_k <- nrow(perf_summ)
         if(!is.null(max.rel.log.mse.increase)){
            warning('Could not find optimal rank (max.rel.log.mse.increase threshold not satisfied).\nReturning last k in k.range')
         }
      } else {
         optimum_k <- perf_summ[(nrow(perf_summ)-1),'k']
      }
   }

   if(verbose>0){ message('Running NMF with optimum rank (=', optimum_k,')...') }
   if(run.with.optimum.k){
      out <- main(k=optimum_k, repeat.num=1, return.perf=F, max.samples=Inf)
   } else {
      out <- list()
   }

   out$perf <- perf
   out$perf_summ <- perf_summ

   return(out)
}


#' Decompose signatures separately for each class
#'
#' @param A A numeric matrix (rows=samples, columns=features)
#' @param response A character or factor vector indicating the class of each sample (i.e. row) in `A`
#' @param nmf.args A list of arguments that can be passed to `nmf()`
#' @param max.sig.dist Signatures will be hierarchical clustered based on pearson correlation
#' distance. The tree will be cut at height==`max.sig.dist`. Signature clusters below this height
#' (i.e. similar signatures) are considered to be redundant. From each cluster, the signature
#' decomposed from the class with the most samples will be kept, and the rest are discarded
#' @param multi.core If TRUE, will use multiple cores
#' @param tmp.dir A path a temporary directory, which if provided, enables resume capability
#' @param verbose Show progress messages? Can be 0,1,2
#'
#' @return A list containing the signature profile matrix (columns=feature, rows=signature) and MSE
#' values from the rank search
#' @export
#'
nmfPerClass <- function(
   A, response, nmf.args=list(), max.sig.dist=0.1,
   multi.core=F, tmp.dir=NULL, verbose=T
){
   if(F){
      A=m_rmd[,-1]
      response=m_rmd$response
      nmf.args=list()
      max.sig.dist=0.1
      multi.core=F
      verbose=2
      tmp.dir='/Users/lnguyen/Desktop/tmp/'
   }

   if(nrow(A) != length(response)){
      stop('No. of rows in `A` must be the same length as `response`')
   }

   A_split <- split(as.data.frame(A), response)

   ## Run NMF for each class --------------------------------
   main <- function(i){
      if(verbose){ message('\nRunning NMF for class: ',i) }

      tmp_out_path <- paste0(tmp.dir,'/nmf_',i,'.rds')
      if(file.exists(tmp_out_path)){
         if(verbose){ message('Reading tmp output: ',tmp_out_path) }
         return(readRDS(tmp_out_path))
      }

      out <- do.call(
         nmf,
         c(list(A_split[[i]]), nmf.args, verbose=(verbose==2))
      )
      saveRDS(out, tmp_out_path)
      return(out)
   }

   if(!multi.core){
      l_nmf_out <- lapply(names(A_split), main)
   }

   else {
      require(foreach)
      require(parallel)
      require(doParallel)

      cores <- parallel::detectCores()
      if(verbose){ message('Multicore: ',cores,' cores') }

      doParallel::registerDoParallel(cores=cores)
      l_nmf_out <- foreach(i=names(A_split)) %dopar% { main(i) }
   }

   names(l_nmf_out) <- names(A_split)

   ## --------------------------------
   if(verbose){ message('Removing redundant signatures...') }
   profiles <- do.call(rbind, lapply(names(l_nmf_out), function(i){
      #i='Anus'
      m <- l_nmf_out[[i]]$H
      rownames(m) <- paste0(i,'.',1:nrow(m))
      return(m)
   }))
   profiles <- profiles/rowSums(profiles)
   profiles[is.na(profiles)] <- 0

   cor_profiles <- cor(t(profiles))

   hc <- hclust(as.dist(1-cor_profiles))
   # plot(hc)
   # abline(h=0)
   # df <- reshape2::melt(cor_profiles)
   # colnames(df) <- c('x','y','sim')
   # df <- df[order(df$sim, decreasing=T),]
   # subset(df, sim<1)

   profile_groups <- split(
      hc$labels,
      cutree(hc, h=max.sig.dist)
   )
   class_freqs <- sort(table(response), decreasing=T)

   profile_whitelist <- lapply(profile_groups, function(i){
      #i=profile_groups[[1]]
      if(length(i)==1){ return(i) }
      sig_class <- gsub('[.]\\d+$','',i)
      class_freqs_ss <- class_freqs[sig_class]
      i[order(class_freqs_ss, decreasing=T)][1]
   })
   profile_whitelist <- unlist(profile_whitelist, use.names=F)

   profiles <- profiles[profile_whitelist,]

   ## --------------------------------
   if(verbose){ message('Summarizing rank search performance...') }
   perf <- do.call(rbind, lapply(names(l_nmf_out), function(i){
      df <- l_nmf_out[[i]]$perf
      df <- cbind(class=i, df)
      return(df)
   }))
   perf$log10_mse_imputed <- log10(perf$mse_imputed)
   perf$mse_imputed <- NULL

   perf_summ <- do.call(rbind, lapply(names(l_nmf_out), function(i){
      df <- l_nmf_out[[i]]$perf_summ
      df <- cbind(class=i, df)
      return(df)
   }))

   ## --------------------------------
   if(verbose){ message('Returning output...') }
   list(
      sig_profiles=profiles,
      perf=perf,
      perf_summ=perf_summ
   )
}











