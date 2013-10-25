#' Calculate Functional Reduncancy
#' 
#' work in progress to calculate FR ...
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param A matrix containing the abundances of the species in x (or presence-absence,
#'  i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#'  species (columns) in a must match the number of species (rows) in x. In addition, 
#'  the species labels in a and x must be identical and in the same order.\code{}
#' @param it number of iterations
#' 
#' @return site vector of sites
#' @return rob_Fpg vector of robustness (indicates reduncancy) using petchey and Gaston index
#' @return rob_FRich vector of robustness using distance-Based indexes (index) Laliberte
#' 
#' @export


FRedundancy <- function(S, A, it = 10){
  rob <- data.frame(site = rownames(A), rob_FDpg = rep(NA,nrow(A)), 
                    rob_FRich = rep(NA,nrow(A)))
  #Note that we remove all species except for two, because you can not 
  #calculate the FD on one species! #I AM PRETTY SURE I FIXED THIS; SO WE CAN UPDATE THE CODE TO REMOVE ALL
  #That may be problematic for AUC's s is a given site
  for (s in 1:nrow(A)){
    #which species are not 0 on the community
    Ano0 <- A[s,-which(A[s,] == 0)] 
    #matrix to store temporal results
    l_funFRich <- matrix(ncol = length(Ano0), nrow = it)
    l_funFD_PG <- matrix(ncol = length(Ano0), nrow = it)
    #create a community matrix with communities and extinct species
    #k is the number of iterations
    for (k in 1:it){ 
      #which species are not 0
      Ano0 <- A[s,-which(A[s,] == 0)]
      Ar <- matrix(nrow = length(Ano0), ncol = nrow(S), 
                   dimnames= list(paste("d", seq(0,length(Ano0)-1,1), sep =""),rownames(S)))
      #create the full community (deleted0)
      Ar[1,] <- as.numeric(A[s,])
      #remove a spe (= 0) one at a time
      for (i in seq(2,length(Ano0)-1,1)){
        #create second community
        Ar[i,] <- Ar[i-1,]
        #remove one random species
        Ar[i,which(colnames(Ar) == names(sample(Ano0,1)))] <- 0
        #remove 0 from the list of species to delete
        Ano0 <- Ar[i,-which(Ar[i,] == 0)]
      } #end i loop
      Ar[i+1,] <- rep(1,ncol(Ar)) #is that correct?? Add all original communities? 
      #calculate FDi for the new communities
      fun <- FDi(S = S, A = Ar, Distance.method= "gower", ord = "podani",
                 w = w, stand.FRic = TRUE, corr = "cailliez" )
      fun[i,8] <- 0
      l_funFRich[k,] <- fun[,8]
      l_funFD_PG[k,] <- fun[,4]
    } #end k loop
    #make the mean
    MeanFunFRich <-  colMeans(l_funFRich)[-(i+1)]
    MeanFunFD_PG <-  colMeans(l_funFD_PG)[-(i+1)]
    #check the variance
    #colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
    #                    twopass=FALSE) {
    #  if (SumSquares) return(colSums(x^2, na.rm, dims))
    #  N <- colSums(!is.na(x), FALSE, dims)
    #  Nm1 <- if (unbiased) N-1 else N
    #  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    #    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
    #  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
    #}
    #colVars(l_funFRich)/sqrt(10)
    #colVars(l_funFD_PG)/sqrt(10) #eyeballing looks like that means are not cancelling very different trajectories
    #calculate the difference
    functFRich <- c()
    functFD_PG <- c()
    for(a in seq(from= nrow(fun)-1, to = 2, by = -1)){
      functFRich[which(a == seq(from= nrow(fun)-1, to = 2,by = -1))] <- 
        abs(MeanFunFRich[a]-MeanFunFRich[a-1])
      functFD_PG[which(a == seq(from= nrow(fun)-1, to = 2,by = -1))] <-
        abs(MeanFunFD_PG[a]-MeanFunFD_PG[a-1])
    } #end a loop
    functFRich[i] <- 0
    functFD_PG[i] <- 0
    #put in bipartite class
    ex.FD <- data.frame(no = c(1:(nrow(fun)-1)), ext.lower = rep(1,nrow(fun)-1), 
                        funct = functFD_PG)
    ex.FD <- as.matrix(ex.FD)
    class(ex.FD) <- "bipartite"
    attr(ex.FD, "exterminated") <- "lower"
    
    ex.FRich <- data.frame(no = c(1:(nrow(fun)-1)), ext.lower = rep(1,nrow(fun)-1), 
                           funct = functFRich) 
    ex.FRich <- as.matrix(ex.FRich)
    class(ex.FRich) <- "bipartite"
    attr(ex.FRich, "exterminated") <- "lower"
    
    #calculate robustness
    rob[s,2] <- robustness(ex.FD)
    #slope.bipartite(ex.FD)
    rob[s,3] <-robustness(ex.FRich)
    #slope.bipartite(ex.FRich)
  } #end s loop
  rob
}
