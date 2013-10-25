#' Dendrogram-Based Functional Diversity weighted Indices used in Clarck et al 2012
#' 
#' Calculate functional trait diversity for a set of communities using Ptchey and Gaston 2002 index 
#' weighted version used in Clarck et al 2012
#' 
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param A matrix containing the abundances of the species in x (or presence-absence,
#'  i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#'  species (columns) in a must match the number of species (rows) in x. In addition, 
#'  the species labels in a and x must be identical and in the same order.\code{}
#' @param w vector listing the weights for the traits in x. Can be missing, 
#' in which case all traits have equal weights.\code{}
#' @param Distance.method metric to calculate the species distance matrix. Only Gower is
#'  implemented. \code{}
#' @param ord character string specifying the method to be used for ordinal traits 
#' (i.e. ordered). "podani" refers to Eqs. 2a-b of Podani (1999), while "metric" 
#' refers to his Eq. 3. Can be abbreviated. See gowdis for more details.\code{}
#' @param Cluster.method Distance method used to produce the tree. UPGMA="average" is 
#' usually giving th ebest results (Podani et al. 2011)\code{}
#' @param stand.x ogical; if all traits are numeric, should they be standardized 
#' to mean 0 and unit variance? If not all traits are numeric, Gower's (1971) 
#' standardization by the range is automatically used; see gowdis for more details.\code{}
#' @param stand.FD logical; should FD be standardized by the max FD, so that FD 
#' is constrained between 0 and 1?\code{}
#' @param Weigthedby character string indicating if should be weighted by `abundance`
#' or `biomassValue`. If biomassValue is in length units for Carabids or bees, 
#' use options `biomasCarabids` or `biomasBees`.\code{}
#' @param  biomassValue numerical vector with body weigh (or length) values for each species
#'  in the same order as species are provided.  \code{}
#'
#' @return comm vector with the name of the community
#' @return n_sp vector listing the number of species for each community
#' @return n_tr vector listing the number of traits used
#' @return qual.FD vector repeating the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' @return FDabund See Clark description.
#' @return FDjointabund See Clark description.
#' 
#' @note This indexes are highly correlated with FDw and FDwcomm, but 
#'  1) the use of recalculated dendograms for each community is advised
#'  2) can not be applied to categorical traits (ignored in the function) 
#' 
#' @export
#' 
#' @examples 
#' ex1 <- FD_Clark(dummy$trait, dummy$abun, Cluster.method = "average", ord = "podani"
#'                     Weigthedby = "abundance")
#' ex1



FD_Clark <- function(S, A, w = NA, Distance.method = "gower", ord= c("podani", "metric"),
                      Cluster.method = c(ward="ward",single="single",complete="complete",
                                         UPGMA="average",UPGMC="centroid",WPGMC="median",
                                         WPGMA="mcquitty"), stand.x = TRUE, stand.FD = FALSE,
                      Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"),
                      biomassValue = NA){
  require(FD)
  require(cluster)
  require(vegan)
  Out <- data.frame(comm = rep(NA,nrow(A)),
                    n_sp = rep(NA,nrow(A)),
                    n_tr = rep(NA,nrow(A)),
                    qual.FD = rep(NA,nrow(A)),
                    FDabund = rep(NA,nrow(A)),
                    FDjointabund = rep(NA,nrow(A))
  )
  Out$comm <- rownames(A)
  Out$n_tr <- ncol(S)
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
  c_distance <- cophenetic(tree)
  Out[, 4] <- rep(cor(D , c_distance), nrow(Out))
  #if Weigthedby is not abundance, transform weight to biomass
  if(Weigthedby != "abundance"){
    if(Weigthedby == "biomasCarabids"){
      biomassValue2 <- Jelaska(biomassValue)
    }
    if(Weigthedby == "biomasBees"){
      biomassValue2 <- Cane(biomassValue)
    }else{
      biomassValue2 <- biomassValue 
    }
    AA <- A*biomassValue2
  }else
    AA <- A
  #create an AFw matrix of relative abundances (/by max)
  AFw <- AA
  for(k in 1:nrow(AA)){
    AFw[k,] <- AA[k,]/max(AA[k,])
  }
  #Calculate FDi for each community
  for(i in 1:nrow(A)){
    
    #Species richness
    Out[i,2] <- length(A[i,-which(A[i,] == 0)]) 
    
    #Calculate FDabund
    #here we have to modify the trait matrix S for each comm!
    if (Distance.method != "gower" & stand.x == TRUE){
      S2 <- scale(S, center = TRUE, scale = TRUE)
    }else{S2 <- S}
    #Select S2 numeric
    num <- sapply(S2, is.numeric)
    Snum <- S2[,num]
    Sabund <- Snum
    if (stand.x == TRUE){
      Sabund <- scale(Sabund, center = TRUE, scale = TRUE)
    }
    #multiply all traits by species relative abundance
    if(is.data.frame(Snum) == TRUE){
      for(k in 1:ncol(Snum)){
        Sabund[,k] <- Snum[,k]*as.numeric(AFw[i,])
      }
    } else{
      for(k in 1:length(Snum)){
        Sabund <- Snum[k]*as.numeric(AFw[i,])
      }
    }
    #Add back categorical variables
    Sabund <- cbind(Sabund, S2[,!num])
    #Remove species not present in the i comm and tweak weigths
    Sabund <- Sabund[which(A[i,] > 0),]
    wabund <- rep(1,ncol(Sabund))
    #Obtain the distance matrix
    if(Distance.method == "gower"){
      #this is fixing an odd error in gowdis, which don't accept ordered factors as numeric if only contains 0 and 1
      is.bin <- function(x) all(x[!is.na(x)] %in% c(0, 1))
      bin.var <- rep(NA, dim(Sabund)[2])
      names(bin.var) <- colnames(Sabund)
      for (j in 1:dim(Sabund)[2]) bin.var[j] <- is.bin(Sabund[, j])
      if(any(bin.var == TRUE)){
        Sabund[,which(bin.var == TRUE)] <- sapply(Sabund[, which(bin.var == TRUE)], as.numeric)
      } else{     #end of the fix
        Dabund <- gowdis(Sabund, w = wabund, ord= ord)
      }
    } else{
      Dabund <- dist(Sabund, method = Distance.method)
    }
    #Obtain the comm based dendrogram
    if(attr(Dabund, "Size") <= 2){
      print("FDabund not calculated for <= 2 species. NA inserted")
      Out[i,5] <- NA
    }else{
      treeabund <- hclust(Dabund, method = Cluster.method)
      #Get the total branch length
      xtreeabund <- Xtree(treeabund)
      #And finnally calculate FDabund
      Acomm <- A[i,which(A[i,] > 0)]
      species_in_Cabund <- ifelse(Acomm>0, 1, 0)
      i.primeCabund <- ifelse(colSums(xtreeabund$H1[which(species_in_Cabund > 0),])>0, 1, 0)
      Out[i,5] <- sum(i.primeCabund*xtreeabund$h2.prime)
    }
    
    ##Calculate FDjointabund
    #modify the D matrix by multiplying by abundance of Species j and Abundance of species a
    Dw <- D
    D2 <- as.matrix(D)
    Dw2 <- as.matrix(Dw)
    for (a in 1:ncol(AFw)){
      for (j in 1:ncol(AFw)){
        Dw2[a,j] <- D2[a,j]*(AFw[i,a]*AFw[i,j]) 
      }
    }
    #remove species not in the comm
    sp_present <- which(A[i,] > 0)
    Dw2 <- Dw2[sp_present, sp_present]
    Dw <- as.dist(Dw2)
    #Obtain the dendrogram for the modified distance matrix
    if(attr(Dw, "Size") <= 2){
      print("FDjointabund not calculated for <= 2 species. NA inserted")
      Out[i,6] <- NA
    }else{
      tree_jointabund <- hclust(Dw, method = Cluster.method)
      #Get the total branch length
      xtree_jointabund <- Xtree(tree_jointabund)
      #calculate the jointabund version
      species_in_Cjointabund <- ifelse(Acomm>0, 1, 0)
      i.primeC_jointabund <- ifelse(colSums(xtree_jointabund$H1[which(species_in_Cjointabund > 0),])>0, 1, 0)
      Out[i,6] <- sum(i.primeC_jointabund*xtree_jointabund$h2.prime)    
    }
  }  
  Out
}
