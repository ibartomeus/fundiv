#' Taxonimic betadiversity
#' 
#' Calculate taxonomic trait beta diversity for two communities using Carvalho et al. 2012 decomposition.
#' 
#' @param c1 vector containing the species names in the first community. NA not tolerated. \code{}
#' @param c2 vector containing the species names in the first community. NA not tolerated. \code{}
#'
#' @return Btot Total Beta diversity
#' @return B_3 Beta diversity due to replacement
#' @return Brich Beta diversity due to richness diferences
#' @return quality the ouptput will print the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- betadiv(c1 = c("sp3", "sp2", "sp1", "sp4", "sp5"), c2 = c("sp6", "sp7", "sp8", "sp4", "sp5"))
#' ex1

beta  <-  function(c1,c2){
  B15 <- function(pm) with(pm,{(b+c)/(a+b+c)})
  B_3  <-  function(pm) with(pm,{2*min(b,c)/(a+b+c)}) 
  Brich  <-  function(pm) with(pm,{abs(b-c)/(a+b+c)})
  pmb  <-  function(A,B) list(b=sum(!(A %in% B)),c=sum(!(B %in% A)),a=sum(B %in% A))
  beta_  <-  B15(pmb(c1,c2))
  beta_3  <-  B_3(pmb(c1,c2))
  beta_rich  <-  Brich(pmb(c1,c2))
  return(list(Btot = beta_, B_3 = beta_3, Brich = beta_rich))
}

