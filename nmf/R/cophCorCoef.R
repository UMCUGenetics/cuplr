#' Calculate the cophenetic correlation coefficient
#'
#' @description The cophenetic measures the stability of the NMF clusters (i.e. signatures).
#'
#' The procedure is below (based on the Brunet NMF paper: https://www.pnas.org/content/101/12/4164):
#' (i) Within each rank, and each iteration, assign each sample a signature/cluster number (=the max
#' contribution signature from the W matrix).
#' (ii) Convert this vector of cluster numbers into a 0/1 matrix of n_samples x n_samples, with a 1
#' indicating 2 samples belong to the sample cluster. This is the 'connectivity matrix' (C in the
#' Brunet paper)
#' (iii) With e.g. 100 NMF repeats, there are 100 connectivity matrices. Take the mean of the
#' connectivity matrices. This is the 'consensus matrix' (C̄ in the paper)
#' (iv) Take 1 - consensus matrix to get the distance matrix, then hierachical cluster on this and
#' calculate the cophenetic matrix with `stats::cophenetic()`
#' (v) The cophenetic coefficient is the pearson correlation of the distance matrix with the
#' cophenetic matrix
#'
#' The above procedure is repeated for each rank.
#'
#' @param clusters A matrix were columns represent samples and rows represent ranks AND iterations.
#' Rownames must be in the form k_{rank number}.rep_{repeat number}
#' @param verbose Show progress messages?
#'
#' @return A numeric vector of cophenetic correlation coefficients
#' @export
#'
cophCorCoef <- function(clusters, verbose=F){

   if(F){
      clusters <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/signatures/04_de_novo/nmf_breast_clean_indel/clusters.txt.gz'
      clusters <- '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/features/signatures/04_de_novo/nmf_breast_clean/clusters.txt.gz'
      verbose=T
   }

   ## --------------------------------
   if(is.vector(clusters) & is.character(clusters)){
      if(verbose){ message('Reading clusters from txt file...') }
      clusters <- read.delim(clusters,header=F, row.names=1)
   }

   clusters <- as.matrix(clusters)

   ## --------------------------------
   k_values <- sapply(strsplit(rownames(clusters),'[.]'),`[[`,1)
   k_values <- as.integer(gsub('^.+_','',k_values))

   k_groups <- split(
      1:nrow(clusters),
      factor(k_values, unique(k_values))
   )

   colnames(clusters) <- NULL
   rownames(clusters) <- NULL

   ## --------------------------------
   if(verbose){ message('Calculating mean connectivity...') }
   calcMeanConnectivity <- function(m){
      #m <- clusters[k_groups[[5]],,drop=F]

      uniq_clusters <- sort(unique(as.vector(m)))
      n_samples <- ncol(m)
      n_iters <- nrow(m)

      l_m_connectivity <- lapply(1:n_iters, function(i){
         #i=1
         cluster_nums <- m[i,]
         m_connectivity <- matrix(0, nrow=n_samples, ncol=n_samples)

         for(j in uniq_clusters){
            #j=5
            in_cluster <- cluster_nums==j
            m_connectivity[in_cluster,in_cluster] <- 1
         }
         return(m_connectivity)
      })

      Reduce("+", l_m_connectivity) / length(l_m_connectivity)
   }

   l_m_mean_connectivity <- lapply(k_groups, function(i){
      #i=k_groups[[5]]
      clusters_ss <- clusters[i,]
      calcMeanConnectivity(clusters_ss)
   })

   ## --------------------------------
   if(verbose){ message('Calculating cophenetic correlation...') }
   coph_cor <- sapply(l_m_mean_connectivity, function(i){
      #i=l_m_mean_connectivity[[5]]
      m_dist <- as.dist(1-i)
      hc <- hclust(m_dist)
      coph <- cophenetic(hc)
      cor(m_dist, coph, method='pearson')
   })

   data.frame(
      k=as.integer(names(coph_cor)),
      coph_cor=coph_cor,
      row.names=NULL
   )

   # library(ggplot2)
   # ggplot(out, aes(x=k, y=coph_cor)) +
   #    geom_point() +
   #    geom_line() +
   #    ylim(0,NA)
}

