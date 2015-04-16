#' Taxonimic betadiversity among several communities
#' 
#' Calculate taxonomic trait beta diversity among all pairwise communities using Carvalho et al. 2012 decomposition
#' and creates a dist object.
#' 
#' @param W matrix containing the abundances of the species per site (or presence-absence,
#' i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. \code{}
#' 
#' @return Btot Total Beta diversity dist object
#' @return B_3 Beta diversity due to replacement dist object
#' @return Brich Beta diversity due to richness diferences dist object
#' @return quality the ouptput will print the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- beta_dist(W = dummy$abun)
#' ex1
#' 

beta_dist = function(W){
  dWN = matrix(NA, ncol = nrow(W), nrow= nrow(W))
  colnames(dWN) = rownames(W)
  rownames(dWN) = rownames(W)
  dBtot = dWN
  dB_3 = dWN
  dBrich = dWN
  for(i in 1:nrow(W)){
    for(j in 1:nrow(W)){
      partition = betadiv(names(which(W[i,] > 0)),names(which(W[j,] > 0)))
      dBtot[i,j] = partition$Btot
      dB_3[i,j] = partition$B_3
      dBrich[i,j] = partition$Brich
    }
  }
  distances = list(Btot = dBtot, B_3 = dB_3, Brich = dBrich)
  distances = lapply(distances,as.dist)
  return(distances)
}


