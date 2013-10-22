FDindexes <- function(S, A, w = NA, Distance.method = "gower", ord= c("podani", "metric"),
                Cluster.method = c(ward="ward",single="single",complete="complete",
                                   UPGMA="average",UPGMC="centroid",WPGMC="median",
                                   WPGMA="mcquitty"),
                stand.x = TRUE, corr = c("sqrt", "cailliez", "lingoes", "none"),
                stand.FRic = FALSE, m = "max", stand.FD = FALSE, scale.RaoQ = FALSE,
                Weigthedby = c("biomasCarabids", "biomasBees", "biomassValue"),
                biomassValue = rep(1,nrow(S))){
  require(FD)
  require(cluster)
  require(vegan)
  Out <- data.frame(comm = rep(NA,nrow(A)),
                    n_sp = rep(NA,nrow(A)),
                    n_tr = rep(NA,nrow(A)),
                    FDpg = rep(NA,nrow(A)),
                    FDw = rep(NA,nrow(A)),
                    FDwcomm = rep(NA,nrow(A)),
                    qual.FD = rep(NA,nrow(A)),
                    FDw_bm = rep(NA,nrow(A)),
                    FDwcomm_bm = rep(NA,nrow(A)),
                    sing.sp = rep(NA,nrow(A)),
                    qual.FRic = rep(NA,nrow(A)),
                    Frich = rep(NA,nrow(A)),
                    Fdis = rep(NA,nrow(A)),
                    Fdis_bm = rep(NA,nrow(A)),
                    Feve = rep(NA,nrow(A)),
                    Fdiv = rep(NA,nrow(A)),
                    RaoQ = rep(NA,nrow(A)),
                    Seve = rep(NA,nrow(A)),
                    Shannon = rep(NA,nrow(A)),
                    Abund = rep(NA,nrow(A)),
                    TotalBiomass = rep(NA,nrow(A)),
                    EvenessBiomass = rep(NA,nrow(A)))
  #calculate FD for metrics using dendrograms
  Outdendro <- FD_dendro(S, A, w = w, Distance.method = Distance.method, ord= ord,
                         Cluster.method = Cluster.method, stand.x = stand.x, stand.FD = stand.FD,
                         Weigthedby = "abundance", biomassValue = NA) 
  Out[,1] <- Outdendro$comm
  Out[,2] <- Outdendro$n_sp
  Out[,3] <- Outdendro$n_tr
  Out[,4] <- Outdendro$FDpg
  Out[,5] <- Outdendro$FDw
  Out[,6] <- Outdendro$FDwcomm
  Out[,7] <- Outdendro$qual.FD
  Outdendro_bm <- FD_dendro(S, A, w = w, Distance.method = Distance.method, ord= ord,
                         Cluster.method = Cluster.method, stand.x = stand.x, stand.FD = stand.FD,
                         Weigthedby = Weigthedby, biomassValue = biomassValue) 
  Out[,8] <- Outdendro_bm$FDw
  Out[,9] <- Outdendro_bm$FDwcomm
  
  #calculate FRichness according to Manson, villager, laliberte et al. using a tweaked version of FD package
  FD <- dbFD_w(S, A, calc.CWM = FALSE, corr = corr, stand.FRic = stand.FRic, 
             ord = ord, m = m, stand.x = stand.x, scale.RaoQ = scale.RaoQ, 
               Weigthedby = "abundance", 
               biomassValue = NA)
  Out[,10] <- FD$sing.sp
  Out[,11] <- rep(FD$qual.FRic, nrow(A))
  Out[,12] <- FD$FRic
  Out[,13] <- FD$FDis
  FD_bm <- dbFD_w(S, A, calc.CWM = FALSE, corr = corr, stand.FRic = stand.FRic, 
               ord = ord, m = m, stand.x = stand.x, scale.RaoQ = scale.RaoQ, 
               Weigthedby = Weigthedby, 
               biomassValue = biomassValue)
  Out[,14] <- FD_bm$FDis
  Out[,15] <- FD$FEve
  Out[,16] <- FD$FDiv
  Out[,17] <- FD$RaoQ
  #calculate species eveness Pielous
  eve <- rep(NA,nrow(A))
  for (k in 1:nrow(A)){
    eve[k] <- diversity(A[k,])/log(specnumber(A[k,]))
  }
  Out[,18] <- eve
  #calculate abundance and shannon
  Out[,19] <- diversity(A)
  Out[,20] <- rowSums(A)
  if(Weigthedby == "biomasCarabids"){
    biomassValue2 <- Jelaska(biomassValue)
  }
  if(Weigthedby == "biomasBees"){
    biomassValue2 <- Cane(biomassValue)
  }else{
    biomassValue2 <- biomassValue 
  }
  Out[,21] <- rowSums(A*biomassValue2)
  #eveness of total biomass
  eve_bm <- rep(NA,nrow(A))
  AA <- A*biomassValue2
  for (k in 1:nrow(AA)){
    eve_bm[k] <- diversity(AA[k,])/log(specnumber(AA[k,]))
  }
  Out[,22] <- eve_bm
  Out
}


#check it works
#FDindexes(dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
#          Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue")

