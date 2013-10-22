#null.FD------
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
                    null_sdFrich = null_sdFrich, pFD = round(pFrich,4))  
  out
}