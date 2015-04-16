#' Functional betadiversity among several communities
#' 
#' Calculate functional trait beta diversity for several communities using Petchey and Gaston 2002 index 
#' and Carvalho et al. 2012 decomposition and creates a dist object.
#' 
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated. All species in W should be present. 
#' Species should be ordered alphabetically.\code{}
#' @param W matrix containing the abundances of the species in S (or presence-absence,
#' i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. 
#' In addition, the species labels in W and S must be identical.\code{}
#'
#' @return Btot Total Betadiversity dist object
#' @return B_3 Beta diversity due to replacement dist object
#' @return Brich Beta diversity due to richness diferences dist object
#' @return quality the ouptput will print the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- betaFD_dist(W = dummy$abun, S = dummy$trait)
#' ex1

betaFD_dist = function(W, S){
  dWN = matrix(NA, ncol = nrow(W), nrow= nrow(W))
  colnames(dWN) = rownames(W)
  rownames(dWN) = rownames(W)
  dBtot = dWN
  dB_3 = dWN
  dBrich = dWN
  treeS <- FDtree(S, w = NA, Distance.method = "gower", ord= "podani",
                  Cluster.method = "average") 
  for(i in 1:nrow(W)){
    for(j in 1:nrow(W)){
      c1 <- names(which(W[i,] > 0))
      c2 <- names(which(W[j,] > 0))      
      partition = betaFD(c1, c2, S, Tree = treeS)
      dBtot[i,j] = partition$Btot
      dB_3[i,j] = partition$B_3
      dBrich[i,j] = partition$Brich
    }
  }
  distances = list(Btot = dBtot, B_3 = dB_3, Brich = dBrich)
  distances = lapply(distances,as.dist)
  return(distances)
}

