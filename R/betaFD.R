#' Functional betadiversity
#' 
#' Calculate functional trait beta diversity for two communities using Petchey and Gaston 2002 index 
#' and Carvalho et al. 2012 decomposition.
#' 
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated. All species in c1 and c2 should be present with no extra species. 
#' Species should be ordered alphabetically.\code{}
#' @param c1 vector containing the species names in the first community. NA not tolerated. In addition, 
#' the species labels in c1 and S must be identical. \code{}
#' @param c2 vector containing the species names in the first community. NA not tolerated. In addition, 
#' the species labels in c2 and S must be identical.\code{}
#'
#' @return Btot Total Betadiversity
#' @return B_3 Beta diversity due to replacement
#' @return Brich Beta diversity due to richness diferences
#' @return quality the ouptput will print the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- betaFD(c1 = c("sp3", "sp2", "sp1", "sp4", "sp5"), c2 = c("sp6", "sp7", "sp8", "sp4", "sp5"), S = dummy$trait)
#' ex1
betaFD <- function(c1,c2,S, Tree = NULL){
  B15 <- function(pm) with(pm,{(b+c)/(a+b+c)})
  B_3  <-  function(pm) with(pm,{2*min(b,c)/(a+b+c)}) 
  Brich  <-  function(pm) with(pm,{abs(b-c)/(a+b+c)})
  b  <-  subset(c1, !c1 %in% c2)
  c <- subset(c2, !c2 %in% c1)
  a <- subset(c2, c2 %in% c1)
  ac <- c(a,c)
  ab <- c(a,b)  
  all <- c(a,b,c)
  A_long <- data.frame(comm = c(rep("a", length(a)),
                               rep("b", length(b)),
                               rep("c", length(c)),
                               rep("ab", length(ab)),
                               rep("ac", length(ac)),
                               rep("all", length(all))),
                  species = c(a,b,c,ab,ac,all)) #note a,b or c can be not present if 0
  A <- table(A_long$comm, A_long$species)
  fd <- FD_dendro2(S, A, Tree, Cluster.method = "average", ord = "podani")
  if(!is.na(fd$qual.FD[1])) print(paste("The quality of the dendogram is", round(fd$qual.FD[1], 2)))
  a = fd[which(fd$comm == "a"), "FDpg"]
  if(!length(a) > 0){
    a = 0
    b = fd[which(fd$comm == "all"), "FDpg"] - fd[which(fd$comm == "c"), "FDpg"]
    c = fd[which(fd$comm == "all"), "FDpg"] - fd[which(fd$comm == "b"), "FDpg"]
    a = fd[which(fd$comm == "b"), "FDpg"] + fd[which(fd$comm == "c"), "FDpg"] - fd[which(fd$comm == "all"), "FDpg"]
  } else {
    b_ = fd[which(fd$comm == "ab"), "FDpg"] - a 
    c_ = fd[which(fd$comm == "ac"), "FDpg"] - a
    a_ = b_ + c_ + a - fd[which(fd$comm == "all"), "FDpg"]
    a = a + a_
    b = b_ - a_
    c = c_ - a_
  }
  branches <- list(a, b, c)
  #recode non presences with 0's
  beta_  <-  B15(branches)
  beta_3  <-  B_3(branches)
  beta_rich  <-  Brich(branches)
  return(list(Btot = round(beta_, 4), B_3 = round(beta_3, 4), Brich = round(beta_rich, 4))) 
}



