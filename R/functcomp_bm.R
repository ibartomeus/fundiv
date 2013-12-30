#' Functional Composition from FD package with the option to weith by biomass
#' 
#' functcomp returns the functional composition of a set of communities, as measured by the community-level weighted means of trait values (CWM; e.g. Lavorel et al. 2008).
#' functcomp_w also allows to weight CWM by biomass
#' 
#' @param same as functcomp\code{}
#' @param Weigthedby character string indicating if should be weighted by `abundance`
#' or `biomassValue`. If biomassValue is in length units for Carabids or bees, 
#' use options `biomasCarabids` or `biomasBees`.\code{}
#'  @param  biomassValue numerical vector with body weigh (or length) values for each species
#'  in the same order as species are provided.  \code{}
#'
#' @return a data frame containing the CWM values (weighted by abundance or by biomass) of each trait for each community.\code{}
#'
#' @export
#' 
#' @examples 
#' ex1 <- functcomp(dummy$trait, dummy$abun, Weigthedby = "biomassValue", 
#' biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
#' ex1

functcomp_bm <- function (x, a, CWM.type = c("dom", "all"), bin.num = NULL,
                          Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"), 
                          biomassValue = NA) {
  if (!is.matrix(x) & !is.data.frame(x)) 
    stop("'x' must be a matrix or a data frame.", "\n")
  else x <- data.frame(x)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.", "\n")
  if (is.null(row.names(x))) 
    stop("'x' must have row names.", "\n")
  else x.n <- row.names(x)
  if (is.null(colnames(a))) 
    stop("'a' must have column names.", "\n")
  else a.n <- colnames(a)
  s.x <- dim(x)[1]
  s.a <- dim(a)[2]
  if (s.x != s.a) 
    stop("Different number of species in 'x' and 'a'.", "\n")
  if (any(x.n != a.n)) 
    stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  com <- dim(a)[1]
  t <- dim(x)[2]
  com.names <- row.names(a)
  sp.names <- row.names(x)
  tr.names <- names(x)
  a[which(is.na(a))] <- 0
  CWM.type <- match.arg(CWM.type)
  is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
  bin.var <- rep(NA, t)
  names(bin.var) <- tr.names
  for (i in 1:t) bin.var[i] <- is.bin(x[, i])
  if (!all(bin.var[bin.num])) 
    stop("'bin.num' points to non-binary variables.\n")
  bin.var[bin.num] <- FALSE
  type <- sapply(x, data.class)
  type[type %in% c("numeric", "integer")] <- "C"
  type[type == "ordered"] <- "O"
  type[type == "factor"] <- "N"
  type[bin.var] <- "B"
  ##IB starts edits
  #if Weigthedby is not abundance, transform weight to biomass
  if(Weigthedby != "abundance"){
    if(Weigthedby == "biomasCarabids"){
      biomassValue2 <- Jelaska(biomassValue)
    }
    if(Weigthedby == "biomasBees"){
      biomassValue2 <- Cane(biomassValue)
    }else{
      biomassValue2 <- biomassValue 
    }
    #now multiply abundances (a) for biomass
    for(i in 1:ncol(a)) a[,i] <- a[,i]*biomassValue2[i]
  }
  ##stop messing up things
  sum.a <- apply(a, 1, sum)
  a <- a/sum.a
  a <- t(a)
  a <- data.frame(a)
  temp <- list()
  for (i in 1:t) {
    if (type[i] == "C") {
      vec <- numeric(com)
      for (j in 1:com) vec[j] <- weighted.mean(x[, i], 
                                               a[, j], na.rm = T)
      temp[[i]] <- matrix(vec, com, 1, dimnames = list(com.names, 
                                                       tr.names[i]))
    }
    else {
      x[, i] <- as.factor(x[, i])
      fac <- data.frame()
      which.dom <- rep(NA, com)
      for (k in 1:com) {
        temp2 <- tapply(a[, k], x[, i], sum)
        fac <- rbind(fac, temp2)
        which.dom[k] <- sample(levels(x[, i])[which(fac[k, 
                                                        ] == max(fac[k, ]))], size = 1)
      }
      colnames(fac) <- paste(tr.names[i], "_", levels(x[, 
                                                        i]), sep = "")
      rownames(fac) <- com.names
      which.dom <- data.frame(which.dom)
      colnames(which.dom) <- tr.names[i]
      rownames(which.dom) <- com.names
      if (CWM.type == "dom") 
        temp[[i]] <- which.dom
      if (CWM.type == "all") 
        temp[[i]] <- fac
    }
  }
  temp <- data.frame(temp)
  return(temp)
}