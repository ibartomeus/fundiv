#' Dendrogram-Based Functional Diversity Indices 
#' 
#' Calculate functional trait diversity for a set of communities using Ptchey and Gaston 2002 index allowing
#' for havind as input a previously calculated tree. 
#'  
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param A matrix containing the abundances of the species in S (or presence-absence,
#' i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#' species (columns) in A must match the number of species (rows) in S. In addition, 
#' the species labels in A and S must be identical and in the same order.\code{}
#' @param Tree a trait dendrogram calculated with FDtree fucntion. 
#' When tree is specified S is ignored. \code{}
#' @param w vector listing the weights for the traits in x. Can be missing, 
#' in which case all traits have equal weights.\code{}
#' @param Distance.method metric to calculate the species distance matrix. Only Gower is
#' implemented. \code{}
#' @param ord character string specifying the method to be used for ordinal traits 
#' (i.e. ordered). "podani" refers to Eqs. 2a-b of Podani (1999), while "metric" 
#' refers to his Eq. 3. See gowdis for more details.\code{}
#' @param Cluster.method Distance method used to produce the tree. UPGMA="average" is 
#' usually giving th ebest results (Podani et al. 2011)\code{}
#' @param stand.x ogical; if all traits are numeric, should they be standardized 
#' to mean 0 and unit variance? If not all traits are numeric, Gower's (1971) 
#' standardization by the range is automatically used; see gowdis for more details.\code{}
#' @param stand.FD logical; should FD be standardized by the max FD, so that FD 
#' is constrained between 0 and 1?\code{}
#'
#' @return comm vector with the name of the community
#' @return n_sp vector listing the number of species for each community
#' @return n_tr vector listing the number of traits used
#' @return FDpg vector listing FDpg (petchey and gaston) for each community
#' @return qual.FD vector repeating the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- FD_dendro2(A = dummy$abun, Tree = FDtree(S = dummy$trait, w = NA, 
#'                    Distance.method = "gower", ord = "podani", Cluster.method = "average"))
#' ex1




FD_dendro2 <- function(S = NULL, A, Tree = NULL, w = NA, Distance.method = "gower", ord= c("podani", "metric"),
                      Cluster.method = c(ward="ward",single="single",complete="complete",
                                         UPGMA="average",UPGMC="centroid",WPGMC="median",
                                         WPGMA="mcquitty"), stand.x = TRUE, stand.FD = FALSE){
  require(FD)
  require(cluster)
  require(vegan)
  Out <- data.frame(comm = rep(NA,nrow(A)),
                    n_sp = rep(NA,nrow(A)),
                    n_tr = rep(NA,nrow(A)),
                    FDpg = rep(NA,nrow(A)),
                    qual.FD = rep(NA,nrow(A))
  )
  Out$comm <- rownames(A)
  Out$n_tr <- ncol(S)
  #richness
  Arich <- as.matrix(A)
  Arich[which(Arich > 0)]  <- 1
  Out$n_sp <- rowSums(Arich, na.rm = TRUE) 
  if(is.null(Tree)){
    if(is.na(w)[1]){w <- rep(1,ncol(S))}
    #Obtain the distance matrix
    if(Distance.method == "gower"){
      D <- gowdis(S, w = w, ord= ord)
    }else{
      if (stand.x == TRUE){
        S2 <- scale(S, center = TRUE, scale = TRUE)
        D <- dist(S2, method = Distance.method)
      }else{
        D <- dist(S, method = Distance.method)
      }
    }
    #Obtain the general dendrogram
    tree <- hclust(D, method = Cluster.method)
    plot(tree)
    #Get the total branch length
    xtree <- Xtree(tree)
    #calculate clustering  performance by using correlation between the cophenetic distance
    c_distance <- cor(D,cophenetic(tree))
    Out[, "qual.FD"] <- rep(c_distance, nrow(Out))
  } else{
    xtree <- Xtree(Tree)
  }
  #Calculate FD for each community
  for(i in 1:nrow(A)){
    species_in_C <- ifelse(A[i,]>0, 1, 0)
    select_xtree <- xtree$H1[which(species_in_C > 0),]
    if(is.array(select_xtree) == TRUE){
      i.primeC <- ifelse(colSums(select_xtree)>0, 1, 0)
    } else{
      i.primeC <- select_xtree
    }
    Out[i,4] <- sum(i.primeC*xtree$h2.prime)
  }
  #standardize FD if needed
  if(stand.FD == TRUE){
    Out[,"FDpg"] <- Out[,4]/max(Out[,4])
  }
  Out
}
