#' Calculate evenness indexes including popular ones and renyi profiles.
#' 
#' Calculate eveness indexes for a set of communities using scripts and the vegan package 
#' and returns a single dataframe. After reading a tone, Pielou's J is included for historical reasons. 
#' Evar and InvSimpson are recomended in X and I Riccota and Avena 2003 recomend to report the spectrum of
#' hill numbers.
#' 
#' @param A matrix containing the abundances of the species in S (or presence-absence,
#' i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#' species (columns) in A must match the number of species (rows) in S. In addition, 
#' the species labels in A and S must be identical and in the same order.\code{}
#' 
#' @param scales vector containing the values of a (hill numbers) to be calculated. 
#' The default is 1, giving you E1,0 index proposed by Riccota and Avena 2003. A more 
#' inclusive exploration would be c(0,0.25,0.5,1,2,4,8,Inf)
#' 
#' @return Ea0 returns a data.frame of class renyi with the selected indices. 
#' @return EinvD returns a vector with EinvD. 
#' @return EJ returns a a vector with Pielou's index. 
#' @return Evar returns a vector with Evar. 
#'
#'
#' @export
#' 
#' @examples 
#' Ea0(A = dummy$abun) 
#' Ea0(A = dummy$abun, scales = c(0.25,0.5,1,2,4,8,Inf)) 
#' EinvD(A = dummy$abun)
#' EJ(A = dummy$abun)
#' Evar(A = dummy$abun)
#' #calculate all of them
#' eve1 <- Eve(A = dummy$abun)
#' eve2 <- Eve(A = dummy$abun, scales = c(0.25,0.5,1,2,4,8,Inf))
#' pairs(eve1)
#' pairs(eve2)


Eve <- function(A, scales = c(1)){
  out <- data.frame(EinvD = rep(NA, nrow(A)),
                    EJ = rep(NA, nrow(A)),
                    Evar = rep(NA, nrow(A)),
                    Ea1 = rep(NA, nrow(A)))
  out$EinvD <- EinvD(A)
  out$EJ <- EJ(A)
  out$Evar <- Evar(A)
  out$Ea1 <- as.numeric(Ea0(A))
  if(length(scales) > 1){
    out[,5:(4+length(scales))] <- Ea0(A, scales = scales)
  }
  rownames(out) <- rownames(A)
  out
}
  
#' @export
#' 
Ea0  <- function(A, scales = c(1)){ 
  require(vegan)
  exp(renyi(A, scales = scales, hill = FALSE)) / 
    exp(renyi(A, scales = c(0), hill = FALSE))
  }

#' @export
EinvD <- function(A) diversity(A, "invsimpson")/specnumber(A)

#' @export
EJ <- function(A){
  eve <- rep(NA,nrow(A))
  for (k in 1:nrow(A)){
    if(specnumber(A[k,]) == 1){
      eve[k] <- NA
    }else{
      eve[k] <- diversity(A[k,])/log(specnumber(A[k,]))
    }
  }
  eve  
}

#' @export
#' 
Evar <- function(A){
  v <- rep(NA, nrow(A)) 
  for(k in 1:nrow(A)) {
    a <- rep(NA, ncol(A)) 
    b <- a  
    d <- a
    for(i in 1:ncol(A)) {
      a[i] <- if(A[k,i] != 0) log(A[k,i]) else 0 
    }
    S <- sum(A[k,] > 0)
    for(i in 1:ncol(A)) {
      b[i] <- a[i]/S
    }
    c <- sum(b)
    for(i in 1:ncol(A)) {
      d[i] <- if(A[k,i] !=0) (a[i]-c)^2/S else 0 
    }
    f <- sum(d)
    v[k] <- (1-2/pi*atan(f))   
  }
  v 
}