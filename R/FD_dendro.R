FD_dendro <- function(S, A, w = NA, Distance.method = "gower", ord= c("podani", "metric"),
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
                    FDpg = rep(NA,nrow(A)),
                    FDw = rep(NA,nrow(A)),
                    FDwcomm = rep(NA,nrow(A)),
                    qual.FD = rep(NA,nrow(A))
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
  Out[, 7] <- rep(cor(D , c_distance), nrow(Out))
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
  #and also weigthed by total
  AFcomm <- AA
  for(k in 1:nrow(AA)){
    AFcomm[k,] <- AA[k,]/max(AA)
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
    #Species richness
    Out[i,2] <- length(A[i,-which(A[i,] == 0)])
    
    ##calculate Fw
    #Substitute all branches where a given species is present (=1) by its weigth 
    xtree.weigths <- xtree$H1
    for(k in 1:nrow(S)){
      xtree.weigths[k,] <- ifelse(xtree$H1[k,] > 0, AFw[i,k], 0)
    }
    #Get the total weigthing for each branch in a vector (i.primeW)
    i.primeW <- c(1:ncol(xtree.weigths))
    for(k in 1:ncol(xtree.weigths)){
      if(sum(xtree.weigths[which(xtree.weigths[,k] > 0),k]) != 0){
        i.primeW[k] <- prod(xtree.weigths[which(xtree.weigths[,k] > 0),k], na.rm = TRUE)
      }else{
        i.primeW[k] <- 0
      }
    }
    #FDw is the sum of the product of i.primeW and h2.prime
    Out[i,5] <- sum(i.primeW*xtree$h2.prime)
    
    ##Calculate FDwcom
    #Substitute all branches where a given species is present (=1) by its weigth 
    #now is raw species numbers... here we divide by the max across comunities
    #that takes into account the abundance.
    xtree.weigths <- xtree$H1
    for(k in 1:nrow(S)){
      xtree.weigths[k,] <- ifelse(xtree$H1[k,] > 0, AFcomm[i,k], 0)
    }
    #Get the total weigthing for each branch in a vector (i.primeW)
    i.primeW <- c(1:ncol(xtree.weigths))
    for(k in 1:ncol(xtree.weigths)){
      if(sum(xtree.weigths[which(xtree.weigths[,k] > 0),k]) != 0){
        i.primeW[k] <- prod(xtree.weigths[which(xtree.weigths[,k] > 0),k], na.rm = TRUE)
      }else{
        i.primeW[k] <- 0
      }
    }
    #FD is the sum of the product of i.primeW and h2.prime
    Out[i,6] <- sum(i.primeW*xtree$h2.prime)
  }
  #standardize FD if needed
  if(stand.FD == TRUE){
    Out[,4] <- Out[,4]/max(Out[,4])
  }
  Out
}




