#' Wrapper to Space and Dendrogram Based Functional Diversity Indices 
#' 
#' Calculate functional trait diversity for a set of communities using FD_dendro 
#' and dbFD and returns a single dataframe
#' 
#' @param S matrix or data frame of functional traits. Traits can be numeric, ordered, 
#' or factor. NAs are tolerated.\code{}
#' @param A matrix containing the abundances of the species in S (or presence-absence,
#' i.e. 0 or 1). Rows are sites and species are columns. NA not tolerated. The number of
#' species (columns) in A must match the number of species (rows) in S. In addition, 
#' the species labels in A and S must be identical and in the same order.\code{}
#' @param w vector listing the weights for the traits in x. Can be missing, 
#' in which case all traits have equal weights.\code{}
#' @param Distance.method metric to calculate the species distance matrix. Only Gower is
#' fully implemented. \code{}
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
#' @param character string specifying the correction method to use when the 
#' species-by-species distance matrix cannot be represented in a Euclidean space. 
#' Options are "sqrt", "cailliez", "lingoes", or "none". Can be abbreviated. 
#' Default is "sqrt". See ‘details’ section.
#' @param stand.FRic logical; should FRic be standardized by the ‘global’ FRic
#' that include all species, so that FRic is constrained between 0 and 1? 
#' @param m the number of PCoA axes to keep as ‘traits’ for calculating FRic 
#' (when FRic is measured as the convex hull volume) and FDiv. Options are: 
#' any integer >1, "min" (maximum number of traits that allows the s >= 2^t 
#' condition to be met, where s is the number of species and t the number of 
#' traits), or "max" (maximum number of axes that allows the s > t condition 
#' to be met). See ‘dbFD’ details section.
#' @param scale.RaoQ logical; should Rao's Q be scaled by its maximal value over
#' all frequency distributions? See divc.
#' @param Weigthedby character string indicating weighted by biomass should be done 
#' on `biomassValue` or corrected first for Carabids or bees. \code{}
#' @param  biomassValue numerical vector with body weigh (or length) values for each species
#' in the same order as species are provided. Default is 1, implying no weightening  \code{}
#' 
#'
#' @return comm vector with the name of the community
#' @return n_sp vector listing the number of species for each community
#' @return n_tr vector listing the number of traits used
#' @return FDpg vector listing FDpg (petchey and gaston) for each community
#' @return FDw vector listing FD weighthed by species relative abundances 
#' in each community
#' @return FDwcomm vector listing FD weighthed by species abundances 
#' across all communities
#' @return qual.FD vector repeating the quality of the dendogram representation.
#' clustering  performance is assessed by the correlation with the cophenetic distance
#' @return FDw_bm FDw vector listing FD weighthed by species relative biomass 
#' in each community
#' @return FDwcomm_bm vector listing FD weighthed by species biomass 
#' across all communities
#' @return sing.sp vector listing the number of functionally singular species 
#' in each community. If all species are functionally different, sing.sp will
#'  be identical to nbsp.
#' @return qual.FRic quality of the reduced-space representation required to 
#' compute FRic and FDiv.
#' @return Frich vector listing the FRic of each community 
#' @return Fdis vector listing the Fdis of each community
#' @return Fdis_bm vector listing the Fdis weighted by biomass of each community
#' @return Feve vector listing the Feve of each community
#' @return Feve_bm vector listing the Feve weighted by biomass of each community
#' @return Fdiv vector listing the Fdiv of each community
#' @return Fdiv_bm vector listing the Fdiv weighted by biomass of each community
#' @return RaoQ vector listing the RaoQ of each community
#' @return RaoQ_bm vector listing the RaoQ weighted by biomass of each community
#' @return Seve vector listing the species eveness (pielous J') of each community
#' @return Shannon vector listing the diversitit Shannon index of each community
#' @return Abund vector listing the Total abundance of each community
#' @return TotalBiomass vector listing the total biomass of each community
#' @return EvenessBiomass vector listing the eveness of biomass distribution
#' of each community
#' @return ShannonBiomass vector listing the Shannon of the diversity weighted by biomass
#' of each community. I am not sure this one makes much sense.
#' 
#' @export
#' 
#' @examples 
#' ex1 <- FDindexes(dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
#'                  Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue",
#'                  biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
#' ex1


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
                    Feve_bm = rep(NA,nrow(A)),
                    Fdiv = rep(NA,nrow(A)),
                    Fdiv_bm = rep(NA,nrow(A)),
                    RaoQ = rep(NA,nrow(A)),
                    RaoQ_bm = rep(NA,nrow(A)),
                    Seve = rep(NA,nrow(A)),
                    Shannon = rep(NA,nrow(A)),
                    Abund = rep(NA,nrow(A)),
                    TotalBiomass = rep(NA,nrow(A)),
                    EvenessBiomass = rep(NA,nrow(A)),
                    ShannonBiomass = rep(NA,nrow(A)))
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
  Out[,16] <- FD_bm$FEve
  Out[,17] <- FD$FDiv
  Out[,18] <- FD_bm$FDiv
  Out[,19] <- FD$RaoQ
  Out[,20] <- FD_bm$RaoQ
  #calculate species eveness Pielous
  eve <- rep(NA,nrow(A))
  for (k in 1:nrow(A)){
    if(specnumber(A[k,]) == 1){
      eve[k] <- "NA"
    }else{
    eve[k] <- diversity(A[k,])/log(specnumber(A[k,]))
    }
  }
  Out[,21] <- eve
  #calculate abundance and shannon
  Out[,22] <- diversity(A)
  Out[,23] <- rowSums(A)
  if(Weigthedby == "biomasCarabids"){
    biomassValue2 <- Jelaska(biomassValue)
  }
  if(Weigthedby == "biomasBees"){
    biomassValue2 <- Cane(biomassValue)
  }else{
    biomassValue2 <- biomassValue 
  }
  AA <- A
  for(i in 1:ncol(A)) AA[,i] <- A[,i]*biomassValue2[i]
  #eveness and Shannon of total biomass
  #Total biomass
  Out[,24] <- rowSums(AA)
  #Eveness_bm
  eve_bm <- rep(NA,nrow(AA))
  for (k in 1:nrow(AA)){
    if(specnumber(AA[k,]) == 1){
      eve[k] <- "NA"
    }else{
    eve_bm[k] <- diversity(AA[k,])/log(specnumber(AA[k,]))
    }
  }
  Out[,25] <- eve_bm
  #Shannon_bm
  Out[,26] <- diversity(AA)  
  Out
}
