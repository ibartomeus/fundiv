#' Calculate Null Model of FD values obtained under random species richness for a given community
#' 
#' work in progress to calculate FR ...
#' 
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param A matrix containing the abundances of the species in x (or presence-absence,
#'  i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#'  species (columns) in a must match the number of species (rows) in x. In addition, 
#'  the species labels in a and x must be identical and in the same order.\code{}
#' @param it number of iterations
#' @param w vector listing the weights for the traits in x. Can be missing, 
#' in which case all traits have equal weights.
#' 
#' @return comm vector listing the evaluated communities
#' @return Rich vector with the observed species richness
#' @return FD  vector with the observed FD index
#' @return null_meanFD vector with the expected mean FD value for that richness 
#' level under random conditions
#' @return null_sdFD vector with the S.E. of the FD value for that richness 
#' level under random conditions
#' @return pFD  p-value comparing FD and its expected null
#' @return Frich vector with the observed FD index
#' @return null_meanFrich vector with the expected mean FD value for that richness 
#' level under random conditions
#' @return null_sdFrich vector with the S.E. of the FD value for that richness 
#' level under random conditions
#' @return pFrich p-value comparing FD and its expected null
#' 
#' 
#' @export

null.FD <- function(S, A, it, w = NA){
  require(GLDEX)
  #select richness levels
  rich_lev <- rowSums(ifelse(A > 0, 1, 0))
  rich_loop <- unique(rich_lev)
  A2 <- A
  #loop thorugh rich loops
  for(r in rich_loop){
    null_com <- matrix(ncol = ncol(A), nrow = it, data = 0)
    colnames(null_com) <- colnames(A)
    rownames(null_com) <- rep(paste("null", r , sep="_"), it)
    for(i in 1:it){
      #fill each vector with 0 and real values according to the null model
      #select i column and add the apropiate number of random species r
      null_com[i, c(sample(c(1:ncol(A)), r, replace = FALSE))] <- 1
      #substitute 1's per abundance values
      for(j in 1:ncol(A)){
        null_com[i,j] <- ifelse(null_com[i,j] == 1, 
                                sample(fun.zero.omit(A[,j]), 1), null_com[i,j])
      } #end loop j
    } #end i loop
    #attch null_com to the list of communities
    A2 <- rbind(A2, null_com)
  } #end r loop
  #Calculate FD in all communities
  fun <- FDi(S = S, A = A2, Distance.method= "gower", ord = "podani",
             w = w, stand.FRic = TRUE, stand.FD = TRUE, corr = "cailliez" )
  #calculate the p-value for FD
  pFD <- c()
  null_meanFD <- c()
  null_sdFD <- c()
  for(k in 1:nrow(A)){
    null_meanFD[k] <-  mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4])
    null_sdFD[k] <-  sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4])
    pFD[k] <- pnorm(fun[k,4], mean = mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4]), 
                    sd = sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),4]))
  } # end k loop
  #calculate the p-value for FRich
  pFrich <- c()
  null_meanFrich <- c()
  null_sdFrich <- c()
  for(k in 1:nrow(A)){
    null_meanFrich[k] <-  mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),9])
    null_sdFrich[k] <-  sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),9])
    pFrich[k] <- pnorm(fun[k,8], mean = mean(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),9]), 
                       sd = sd(fun[tail(which(fun$n_sp == fun$n_sp[k]),it),9]))
  } # end k loop
  #output in table format
  out <- data.frame(comm = fun$comm[1:nrow(A)], Rich = fun$n_sp[1:nrow(A)] ,FD = fun$FD_PG[1:nrow(A)], 
                    null_meanFD = null_meanFD, null_sdFD = null_sdFD, pFD = round(pFD,4), 
                    Frich = fun$Frich[1:nrow(A)], null_meanFrich = null_meanFrich,
                    null_sdFrich = null_sdFrich, pFrich = round(pFrich,4))  
  out
}