fundiv: analyzing functional trait diversity
========================================================

This is a wrapper to Petchey and Gaston FD indexes and Laliberte FD package. Hence most code is borrowed from them. It has some additions, most notably the dendogram measures can also be weighted by abundance as implemented in Gagic, Bartomeus et al. (Submitted) paper. There is also a function to calculate several Evennes indexes.

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


You may also want indexes based in trait space in FD package. You can run 'dbFD {FD}',
or the wrapper below which alos includes other diversity metrics.

```{r}
FD_all <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "abundance")
FD_all
```

A nice addition if that you may want the indexes weigthed by biomass, not abundance. For bees and carabids, a function to convert body length to body mass is provided (see ?length_to_mass)

```{r}
FD_all_bm <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue",
                  biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))

FD_all_bm
```

We can see that both families of indexes are correlated

```plot(FD_all$FDpg ~ FD_all$Frich)```

Note that Functional richness indexes are highly correlated with richness, hence you may also want to know if those indexes are higher or lower than expected by its richness levels (See Rader et al. 2014 D&D in press for details)

```{r}
null <- null.FD(S = dummy$trait, A = dummy$abun, it = 100, w = NA)
null
```

No much significant results in this dataset, but is not surprising given:

```plot(FD_all$FDpg ~ FD_all$n_sp)```

You can calculate standardized FD indexes as in Radet et al 2014 by 

```{r}
(null$FD - null$null_meanFD) / null$null_sdFD
```

Indexes implemented in Clark et al. 2012 are also available (`FD_Clark`), but in my opinion they perform very similar that the proposed FDw, but are computational time consuming because they need to 
create a dendogram for each community, instead of using only a general dendogram containing all species. This last approach is better according to the literature...

Lastly I added the calculation of some evenness indexes. See ?Eve for details.

```{r}
eve1 <- Eve(A = dummy$abun)
eve2 <- Eve(A = dummy$abun, scales = c(0.25,0.5,2,4,8,Inf))
pairs(eve1)
pairs(eve2)
```


