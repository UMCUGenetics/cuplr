#' Run NNLM::nnmf() and calculate MSE
#'
#' @rdname nmf
#' @export
runNmf <- function(
   A, k=1, repeat.num=1, max.samples=NULL,
   perf.metrics=c('mse_imputed', 'mse_perm'),
   impute.prop=0.1, return.perf=1,
   seed=1, verbose=F
){

   # if(F){
   #    k=5
   #    repeat.num=1
   #    max.samples=NULL
   #    perf.metrics=c('mse_imputed', 'mse_perm')
   #    impute.prop=0.1
   #    return.perf=1
   #    seed=1
   #    verbose=T
   # }

   ## Init --------------------------------
   i_seed <- as.integer(paste0(seed, repeat.num))
   set.seed(i_seed)

   A <- as.matrix(A)
   if(!is.numeric(A)){ stop('`A` must be a numeric matrix or dataframe') }

   args.nnmf <- list(
      loss='mkl', max.iter=2000L,
      verbose=if(verbose==2){ 2 } else { 0 }
   )

   ## Reduce number of samples to reduce run time
   if(!is.null(max.samples)){
      if(nrow(A)>max.samples){
         A <- A[sample(nrow(A), max.samples),]
      }
   }

   ## Initialize performance output
   perf <- data.frame(
      k=k,
      rep=repeat.num,
      mse=NA,
      mse_imputed=NA,
      mse_perm=NA
   )

   ## --------------------------------
   if('mse_imputed' %in% perf.metrics){

      if(verbose){ message('Calculating MSE on imputed values...') }

      ## Randomly replace some values with NA and do nnmf(), which will impute the NA values
      ## https://mran.microsoft.com/snapshot/2017-01-23/web/packages/NNLM/vignettes/Fast-And-Versatile-NMF.html
      ## See 'Determine rank k via missing value imputation'

      A.with_na <- A
      ind <- sample(length(A.with_na), impute.prop*length(A.with_na))
      A.with_na[ind] <- NA

      nmf_out.with_na <- do.call(
         what=NNLM::nnmf,
         args=c(list(A=A.with_na, k=k), args.nnmf)
      )

      A.with_na.reconstr <- nmf_out.with_na$W %*% nmf_out.with_na$H
      perf$mse_imputed <- mean( (A.with_na.reconstr[ind] - A[ind])^2 )

      rm(A.with_na, A.with_na.reconstr, ind)
   }

   ## --------------------------------
   if('mse_perm' %in% perf.metrics){

      if(verbose){ message('Calculating MSE on permuted data...') }

      ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2623306/
      ## Copy paste from paper:
      ## the residual error of A, RE = |A − WH| is computed for different choices of rank and
      ## compared to the residual error of Aperm, RE = |Aperm − W(Aperm)H(Aperm)|, where Aperm
      ## denotes the expression matrix A with the rows (genes) permuted for every column (sample).
      ## The slope in a plot of the residual error versus rank k is a measure of how much
      ## information is lost as the rank decreases. If the slope of the residual error of A is larger
      ## than that of Aperm this indicates information present in the original data set. We thus
      ## identified the smallest rank value which still preserves additional information compared to
      ## Aperm.

      A.perm <- A
      A.perm <- apply(A.perm, 1, function(i){
         i[ sample.int(length(i), replace=F) ]
      })
      A.perm <- t(A.perm)

      nmf_out.perm <- do.call(
         what=NNLM::nnmf,
         args=c(list(A=A.perm, k=k), args.nnmf)
      )

      A.perm.recontr <- nmf_out.perm$W %*% nmf_out.perm$H
      perf$mse_perm <- mean( (A.perm.recontr - A.perm)^2 )

      rm(A.perm, A.perm.recontr)
   }

   ## --------------------------------
   if(return.perf){ return(perf) }

   ## --------------------------------
   if(verbose){ message('Running NMF on full data...') }
   nmf_out.full <- do.call(
      what=NNLM::nnmf,
      args=c(list(A=A, k=k), args.nnmf)
   )

   ##
   perf$mse <- tail(nmf_out.full$mse, 1)

   ##
   nmf_out.full$perf <- perf
   nmf_out.full$clusters <- max.col(nmf_out.full$W)

   return(nmf_out.full)
}

####################################################################################################
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
nmfRankSearch <- function(
   A,
   k.range=1:10, repeats=10, max.samples=100, impute.prop=0.1, max.rel.log.mse.increase=0.002,
   run.with.optimum.k=T, seed=1, verbose=0
){
   if(F){
      A=as.matrix(contexts$snv)
      k.range=2:3
      repeats=3
      max.samples=100
      impute.prop=0.1
      max.rel.log.mse.increase=0.002
      seed=1
      verbose=T
   }

   ## Init --------------------------------
   set.seed(seed)

   A <- as.matrix(A)
   if(!is.numeric(A)){ stop('`A` must be a numeric matrix or dataframe') }

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

      i_perf <- runNmf(
         A=A, k=k, repeat.num=repeat_num, return.perf=T, perf.metrics='mse_imputed',
         impute.prop=impute.prop, max.samples=max.samples, seed=seed, verbose=verbose
      )
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
      out <- runNmf(
         A=A, k=optimum_k, repeat.num=1, return.perf=F, perf.metrics=NULL,
         max.samples=NULL, impute.prop=impute.prop, seed=seed, verbose=verbose
      )
   } else {
      out <- list()
   }

   out$perf <- perf
   out$perf_summ <- perf_summ

   return(out)
}
