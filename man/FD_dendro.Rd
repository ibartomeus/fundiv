\name{FD_dendro}
\alias{FD_dendro}
\title{Dendrogram-Based Functional Diversity Indices}
\usage{
  FD_dendro(S, A, w = NA, Distance.method = "gower",
    ord = c("podani", "metric"),
    Cluster.method = c(ward = "ward", single = "single", complete = "complete", UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA = "mcquitty"),
    stand.x = TRUE, stand.FD = FALSE,
    Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"),
    biomassValue = NA)
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
  distance matrix. Only Gower is implemented. \code{}}

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

  \item{Weigthedby}{character string indicating if should
  be weighted by `abundance` or `biomassValue`. If
  biomassValue is in length units for Carabids or bees, use
  options `biomasCarabids` or `biomasBees` to autmatically
  convert to mass.\code{}}

  \item{biomassValue}{numerical vector with body weigh (or
  length) values for each species in the same order as
  species are provided.  \code{}}
}
\value{
  comm vector with the name of the community

  n_sp vector listing the number of species for each
  community

  n_tr vector listing the number of traits used

  FDpg vector listing FDpg (petchey and gaston) for each
  community

  FDw vector listing FD weighthed by species relative
  abundances/biomass in each community

  FDwcomm vector listing FD weighthed by species
  abundances/biomass across all communities

  qual.FD vector repeating the quality of the dendogram
  representation. clustering performance is assessed by the
  correlation with the cophenetic distance
}
\description{
  Calculate functional trait diversity for a set of
  communities using Ptchey and Gaston 2002 index and its
  weighted version used in Gagig et al. In prep.
}
\examples{
ex1 <- FD_dendro(S = dummy$trait, A = dummy$abun, Cluster.method = "average", ord = "podani",
                    Weigthedby = "abundance")
ex1
}

