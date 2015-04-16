#' Cretae a dendrogram for use in FD analysis 
#' 
#' Calculate dendrogram and extract branch lengths for use in FD analysis. 
#'  
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param w vector listing the weights for the traits in x. Can be missing, 
#' in which case all traits have equal weights.\code{}
#' @param Distance.method metric to calculate the species distance matrix. Only Gower is
#' implemented. \code{}
#' @param ord character string specifying the method to be used for ordinal traits 
#' (i.e. ordered). "podani" refers to Eqs. 2a-b of Podani (1999), while "metric" 
#' refers to his Eq. 3. See gowdis for more details.\code{}
#' @param Cluster.method Distance method used to produce the tree. UPGMA="average" is 
#' usually giving th ebest results (Podani et al. 2011)\code{}
#'
#' @return an xtree object
#' @return The queality of the tree is printed. The quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- FDtree(S = dummy$trait, w = NA, 
#'                    Distance.method = "gower", ord = "podani", Cluster.method = "average")
#' ex1
FDtree <- function(S, w = NA, Distance.method = "gower", ord= c("podani", "metric"),
                   Cluster.method = c(ward="ward",single="single",complete="complete",
                                      UPGMA="average",UPGMC="centroid",WPGMC="median",
                                      WPGMA="mcquitty")){
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
  print(paste("The quality of the dendogram is", round(c_distance, 2)))
  xtree
}