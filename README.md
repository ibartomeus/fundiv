fundiv: analyzing functional trait diversity
========================================================

This is a wrapper to Petchey and Gaston FD indexes and Laliberte FD package. Hence most code is borrowed from them. It has some additions to allow the dendogram to be weighted by abundance implemented in Gagic, Bartomeus et al. (Submitted) paper.

This version is in development, so please report bugs, etc...

To install the package run:

```{r}
install.packages("devtools")
require(devtools)
install_github("fundiv", "ibartomeus")
require(fundiv)
```

To calculate dendogram indexes:

```{r}
FD_dendro <- FD_dendro(S = dummy$trait, A = dummy$abun, Cluster.method = "average", ord = "podani",
                    Weigthedby = "abundance")
FD_dendro
```

The dendrogram is ploted and returns a dataframe with the indexes.


You may also want indexes based in trait space in FD package.

```{r}
FD_all <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "abundance")
FD_all
```

A nice addition if that you may want the indexes weigthed by biomass, not abundance.

```{r}
FD_all_bm <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue",
                  biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))

FD_all_bm
```

We can see that both families of indexes are correlated

```plot(FD_all$FDpg ~ FD_all$Frich)```

Note that Functional richness indexes are highly correlated with richness, hence you may also want to know if those indexes are higher or lower than expected by its richness levels (See Rader et al. 2014 Submitted for details)

```{r}
null.FD(S = dummy$trait, A = dummy$abun, it = 100, w = NA)
```

No much significant results in this dataset, but is not surprising given:

```plot(FD_all$FDpg ~ FD_all$n_sp)```

Indexes implemented in Clark et al. 2012 are also available (`FD_Clark`), but in my opinion they perform worst that the proposed FDw (i'll expand on that later)

