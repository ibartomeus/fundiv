\name{FDindexes}
\alias{FDindexes}
\title{Wrapper to Space and Dendrogram Based Functional Diversity Indices}
\usage{
FDindexes(S, A, w = NA, Distance.method = "gower", ord = c("podani",
  "metric"), Cluster.method = c(ward = "ward", single = "single", complete =
  "complete", UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA =
  "mcquitty"), calc.FRic = TRUE, stand.x = TRUE, corr = c("sqrt",
  "cailliez", "lingoes", "none"), stand.FRic = FALSE, m = "max",
  stand.FD = FALSE, scale.RaoQ = FALSE, Weigthedby = c("biomasCarabids",
  "biomasBees", "biomassValue"), biomassValue = rep(1, nrow(S)))
}
\arguments{
  \item{S}{matrix or data frame of functional traits.
  Traits can be numeric, ordered, or factor. NAs are
  tolerated.\code{}}

  \item{A}{matrix containing the abundances of the species
  in S (or presence-absence, i.e. 0 or 1). Rows are sites
  and species are columns. NA not tolerated. The number of
  species (columns) in A must match the number of species
  (rows) in S. In addition, the species labels in A and S
  must be identical and in the same order.\code{}}

  \item{w}{vector listing the weights for the traits in x.
  Can be missing, in which case all traits have equal
  weights.\code{}}

  \item{Distance.method}{metric to calculate the species
  distance matrix. Only Gower is fully implemented.
  \code{}}

  \item{ord}{character string specifying the method to be
  used for ordinal traits (i.e. ordered). "podani" refers
  to Eqs. 2a-b of Podani (1999), while "metric" refers to
  his Eq. 3. See gowdis for more details.\code{}}

  \item{Cluster.method}{Distance method used to produce the
  tree. UPGMA="average" is usually giving th ebest results
  (Podani et al. 2011)\code{}}

  \item{stand.x}{ogical; if all traits are numeric, should
  they be standardized to mean 0 and unit variance? If not
  all traits are numeric, Gower's (1971) standardization by
  the range is automatically used; see gowdis for more
  details.\code{}}

  \item{stand.FD}{logical; should FD be standardized by the
  max FD, so that FD is constrained between 0 and
  1?\code{}}

  \item{character}{string specifying the correction method
  to use when the species-by-species distance matrix cannot
  be represented in a Euclidean space. Options are "sqrt",
  "cailliez", "lingoes", or "none". Can be abbreviated.
  Default is "sqrt". See ‘details’ section.}

  \item{stand.FRic}{logical; should FRic be standardized by
  the ‘global’ FRic that include all species, so that FRic
  is constrained between 0 and 1?}

  \item{m}{the number of PCoA axes to keep as ‘traits’ for
  calculating FRic (when FRic is measured as the convex
  hull volume) and FDiv. Options are: any integer >1, "min"
  (maximum number of traits that allows the s >= 2^t
  condition to be met, where s is the number of species and
  t the number of traits), or "max" (maximum number of axes
  that allows the s > t condition to be met). See ‘dbFD’
  details section.}

  \item{scale.RaoQ}{logical; should Rao's Q be scaled by
  its maximal value over all frequency distributions? See
  divc.}

  \item{Weigthedby}{character string indicating weighted by
  biomass should be done on `biomassValue` or corrected
  first for Carabids or bees. \code{}}

  \item{biomassValue}{numerical vector with body weigh (or
  length) values for each species in the same order as
  species are provided. It can also be a matrix or data
  frame with one mass value for each community and species
  (both communities and species arranged like in A).
  Default is 1, implying no weightening \code{}}
}
\value{
comm vector with the name of the community

n_sp vector listing the number of species for each
community

n_tr vector listing the number of traits used

FDpg vector listing FDpg (petchey and gaston) for each
community

FDw vector listing FD weighthed by species relative
abundances in each community

FDwcomm vector listing FD weighthed by species abundances
across all communities

qual.FD vector repeating the quality of the dendogram
representation. clustering performance is assessed by the
correlation with the cophenetic distance

FDw_bm FDw vector listing FD weighthed by species relative
biomass in each community

FDwcomm_bm vector listing FD weighthed by species biomass
across all communities

sing.sp vector listing the number of functionally singular
species in each community. If all species are functionally
different, sing.sp will be identical to nbsp.

qual.FRic quality of the reduced-space representation
required to compute FRic and FDiv.

Frich vector listing the FRic of each community

Fdis vector listing the Fdis of each community

Fdis_bm vector listing the Fdis weighted by biomass of each
community

Feve vector listing the Feve of each community

Feve_bm vector listing the Feve weighted by biomass of each
community

Fdiv vector listing the Fdiv of each community

Fdiv_bm vector listing the Fdiv weighted by biomass of each
community

RaoQ vector listing the RaoQ of each community

RaoQ_bm vector listing the RaoQ weighted by biomass of each
community

Seve vector listing the species eveness (pielous J') of
each community

Shannon vector listing the diversitit Shannon index of each
community

Abund vector listing the Total abundance of each community

TotalBiomass vector listing the total biomass of each
community

EvenessBiomass vector listing the eveness of biomass
distribution of each community

ShannonBiomass vector listing the Shannon of the diversity
weighted by biomass of each community. I am not sure this
one makes much sense.
}
\description{
Calculate functional trait diversity for a set of
communities using FD_dendro and dbFD and returns a single
dataframe
}
\examples{
ex1 <- FDindexes(dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani",
                 Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue",
                 biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
ex1
}

