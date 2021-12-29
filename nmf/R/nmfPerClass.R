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
         nmfRankSearch,
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

   if(F){
      tmp.dir <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/models/0.18a_newLabels/final/tmp/nmf/'
      tmp_out_paths <- list.files(tmp.dir,full.names=T)
      l_nmf_out <- lapply(tmp_out_paths, readRDS)
      names(l_nmf_out) <- sub(
         '^nmf_','',
         sub('[.]rds','',basename(tmp_out_paths))
      )
   }


   ## --------------------------------
   if(verbose){ message('Removing redundant signatures...') }
   ## Name signatures like so: {cancer_type}.{denovo_sig_number}
   profiles <- do.call(rbind, lapply(names(l_nmf_out), function(i){
      #i='Anus'
      m <- l_nmf_out[[i]]$H
      rownames(m) <- paste0(i,'.',1:nrow(m))
      return(m)
   }))

   ## Scale signature profiles to sum to 1. (The 'H' matrix from NMF is not by default properly scaled)
   profiles <- profiles/rowSums(profiles)
   profiles[is.na(profiles)] <- 0

   ## Clustering based on pearson correlation
   cor_profiles <- cor(t(profiles))
   distances <- as.dist(1-cor_profiles)
   hc <- hclust(distances)
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

   ## For similar signature profiles, select the one from largest cancer type cohort.
   ## The logic is that the signature extraction from a larger cohort gives a more stable signature.
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

   # if(F){
   #    #n_clust_range <- 2:min(nrow(cor_profiles),100)
   #    n_clust_range <- 2:46
   #    avg_sil_width <- sapply(n_clust_range, function(i){
   #       #i=3
   #       clusters <- cutree(hc, k=i)
   #       sil <- cluster::silhouette(clusters, distances)
   #       mean(sil[,3])
   #    })
   #
   #    require(ggplot2)
   #    pd <- data.frame(k=n_clust_range, avg_sil_width)
   #    ggplot(pd, aes(x=k, y=avg_sil_width)) +
   #       geom_point() +
   #       geom_line() +
   #       theme_bw()
   # }

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











