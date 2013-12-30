#' Distance-Based Functional Diversity Indices from FD package with the option to weith by biomass
#' 
#' dbFD implements a flexible distance-based framework to compute multidimensional functional diversity (FD) indices. dbFD returns the three FD indices of Villéger et al. (2008): functional richness (FRic), functional evenness (FEve), and functional divergence (FDiv), as well functional dispersion (FDis; Laliberté and Legendre 2010), Rao's quadratic entropy (Q) (Botta-Dukát 2005), a posteriori functional group richness (FGR) (Petchey and Gaston 2006), and the community-level weighted means of trait values (CWM; e.g. Lavorel et al. 2008). Some of these FD indices consider species abundances. dbFD includes several options for flexibility.
#' dbFD_w also allows to weight FDis by biomass
#' 
#' @param same as dbFD\code{}
#' @param Weigthedby character string indicating if should be weighted by `abundance`
#' or `biomassValue`. If biomassValue is in length units for Carabids or bees, 
#' use options `biomasCarabids` or `biomasBees`.\code{}
#'  @param  biomassValue numerical vector with body weigh (or length) values for each species
#'  in the same order as species are provided.  \code{}
#'
#' @return same as dbFD\code{}
#' @return FDis vector listing the FDis of each community weighted by abundance or biomass
#'
#' @export
#' 
#' @examples 
#' ex1 <- dbFD_w(x = dummy$trait, a = dummy$abun, Weigthedby = "biomassValue", 
#' biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
#' ex1

dbFD_w <- function (x, a, w, w.abun = TRUE, stand.x = TRUE, 
                    ord = c("podani","metric"), asym.bin = NULL, 
                    corr = c("sqrt", "cailliez", "lingoes", "none"), 
                    calc.FRic = TRUE, m = "max", stand.FRic = FALSE, 
                    scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward", 
                    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, 
                    km.crit = c("calinski", "ssi"), calc.CWM = TRUE, 
                    CWM.type = c("dom", "all"), calc.FDiv = TRUE, dist.bin = 2, 
                    print.pco = FALSE, messages = TRUE, 
                    Weigthedby = c("abundance", "biomasCarabids", "biomasBees",
                                   "biomassValue"), biomassValue = NA) {
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  CWM.type <- match.arg(CWM.type)
  km.crit <- match.arg(km.crit)
  if (!is.logical(messages)) 
    stop("'messages' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.FRic)) 
    stop("'stand.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.x)) 
    stop("'stand.x' must be TRUE or FALSE.", "\n")
  if (!is.logical(w.abun)) 
    stop("'w.abun' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FRic)) 
    stop("'calc.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FDiv)) 
    stop("'calc.FDiv' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FGR)) 
    stop("'calc.FGR' musts be TRUE or FALSE.", "\n")
  if (!is.logical(calc.CWM)) 
    stop("'calc.CWM' must be TRUE or FALSE.", "\n")
  if (!is.logical(scale.RaoQ)) 
    stop("'scale.RaoQ' must be TRUE or FALSE.", "\n")
  if (!is.logical(print.pco)) 
    stop("'print.pco' must be TRUE or FALSE.", "\n")
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- dim(x)[1]
    t.x <- dim(x)[2]
    if (is.null(row.names(x))) 
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x))) 
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels"))) 
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (missing(a)) {
    ab.names <- list("Community1", x.rn)
    a <- matrix(1, 1, s.x, dimnames = ab.names)
  }
  else {
    if (is.matrix(a) | is.data.frame(a)) {
      s.a <- dim(a)[2]
      ab.t <- t(a)
      if (is.null(row.names(ab.t))) 
        stop("'a' must have column names.", "\n")
      else ab.t.row <- row.names(ab.t)
      a <- as.matrix(a)
    }
    if (is.vector(a)) {
      s.a <- length(a)
      if (is.null(names(a))) 
        stop("'a' must have names.", "\n")
      else ab.t.row <- names(a)
      ab.names <- list("Community1", ab.t.row)
      a <- matrix(a, 1, s.a, dimnames = ab.names)
    }
    if (s.x != s.a) 
      stop("Different number of species in 'x' and 'a'.", 
           "\n")
    if (any(ab.t.row != x.rn)) 
      stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
           "\n")
  }
  a <- as.matrix(a)
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    stop("At least one community has zero-sum abundances (no species).", 
         "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", 
         "\n")
  if (!missing(w) & is.dist.x) 
    stop("When 'x' is a distance matrix, 'w' should be left missing.", 
         "\n")
  if (!missing(w) & !is.dist.x) {
    if (!is.numeric(w) | length(w) != t.x) 
      stop("'w' should be a numeric vector of length = number of traits.", 
           "\n")
    else w <- w/sum(w)
  }
  if (missing(w)) 
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) 
        x[, x.class == "character"] <- as.factor(x[, 
                                                   x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1])) 
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          x.rn <- x.rn[-pos.NA]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)])) 
            stop("'dist.bin' must be an integer between 1 and 10.", 
                 "\n")
          x.dist <- dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = T, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.", 
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x))) 
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).", 
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist))) 
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0)) 
      stop("At least one species has no trait data.", "\n")
  }
  c <- dim(a)[1]
  if (!w.abun) 
    for (h in 1:c) {
      abpos <- which(a[h, ] > 0)
      a[h, abpos] <- 1
    }
  attr(x.dist, "Labels") <- x.rn
  if (is.euclid(x.dist)) 
    x.dist2 <- x.dist
  if (!is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- lingoes(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", 
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", 
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2)) 
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", 
             "\n")
      if (is.euclid(x.dist2)) 
        if (messages) 
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", 
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages) 
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", 
            "\n")
    }
  }
  x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  traits <- x.pco$li
  nb.sp <- numeric(c)
  for (i in 1:c) {
    sp.pres <- which(a[i, ] > 0)
    traits.sp.pres <- traits[sp.pres, , drop = F]
    traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < 
                     tol] <- 0
    nb.sp[i] <- nrow(unique(traits.sp.pres))
  }
  names(nb.sp) <- row.names(a)
  min.nb.sp <- min(nb.sp)
  if (min.nb.sp < 3) 
    if (messages) 
      cat("FEVe: Could not be calculated for communities with <3 functionally singular species.", 
          "\n")
  if (min.nb.sp < 2) 
    if (messages) 
      cat("FDis: Equals 0 in communities with only one functionally singular species.", 
          "\n")
  if (calc.FRic) {
    x.class2 <- sapply(x, data.class)
    if (all(x.class2 == "factor" | x.class2 == "ordered")) {
      if (length(x.class2) == 1 & x.class2[1] == "ordered") {
        traits.FRic1 <- rank(x[, 1])
        names(traits.FRic1) <- x.rn
        traits.FRic <- data.frame(traits.FRic1)
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one ordinal trait present in 'x'. FRic was measured as the range of the ranks, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when 'x' is a single ordinal trait.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      else {
        traits.FRic <- x
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only categorical and/or ordinal trait(s) present in 'x'. FRic was measured as the number of unique trait combinations, NOT as the convex hull volume.", 
              "\n")
        if (stand.FRic) 
          FRic.all <- nrow((unique(traits.FRic)))
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when only categorical and/or ordinal trait(s) present in 'x'.", 
                "\n")
        }
      }
    }
    else {
      if (x.pco$nf == 1) {
        traits.FRic <- x.pco$li
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      if (x.pco$nf > 1) {
        warning <- FALSE
        m.max <- min.nb.sp - 1
        if (m == "min") {
          warning <- TRUE
          if (min.nb.sp < 4) {
            nb.sp2 <- nb.sp[nb.sp > 3]
            m.min <- floor(log2(min(nb.sp2)))
            if (messages) 
              cat("FRic: To respect s >= 2^t, FRic could not be calculated for communities with <4 functionally singular species.", 
                  "\n")
          }
          else m.min <- floor(log2(min.nb.sp))
        }
        else {
          if (min.nb.sp < 3) {
            nb.sp2 <- nb.sp[nb.sp > 2]
            m.max <- min(nb.sp2) - 1
            if (messages) 
              cat("FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species.", 
                  "\n")
          }
          else m.max <- m.max
        }
        if (is.numeric(m) & m <= 1) 
          stop("When 'm' is an integer, it must be >1.", 
               "\n")
        if (is.numeric(m) & m > m.max) 
          m <- m.max
        if (m == "min") 
          m <- m.min
        if (m == "max") 
          m <- m.max
        if (!is.numeric(m) & m != "min" & m != "max") 
          stop("'m' must be an integer >1, 'min', or 'max'.", 
               "\n")
        if (m < x.pco$nf) {
          traits.FRic <- x.pco$li[, 1:m]
          if (x.pco$nf - m == 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last PCoA axis (out of", 
                  x.pco$nf, "in total) was removed.", "\n")
          if (x.pco$nf - m > 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last", 
                  x.pco$nf - m, "PCoA axes (out of", x.pco$nf, 
                  "in total) were removed.", "\n")
          if (is.euclid(x.dist)) {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr != "none") {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (based on corrected distance matrix) =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr == "none") {
            delta <- -0.5 * bicenter.wt(x.dist * x.dist)
            lambda <- eigen(delta, symmetric = TRUE, 
                            only.values = TRUE)$values
            sum.m <- sum(lambda[1:m])
            sum.n <- sum(lambda)
            lambda.neg <- c(lambda[lambda < 0])
            max.neg <- abs(min(lambda.neg))
            qual.FRic <- (sum.m + (length(lambda[1:m]) * 
                                     max.neg))/(sum.n + ((length(lambda) - 1) * 
                                                           max.neg))
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) =", 
                  qual.FRic, "\n")
          }
        }
        if (m >= x.pco$nf) {
          qual.FRic = 1
          traits.FRic <- x.pco$li
          if (x.pco$nf == 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. The 2 PCoA axes were kept as 'traits'.", 
                  "\n")
          if (x.pco$nf > 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. All", 
                  x.pco$nf, "PCoA axes were kept as 'traits'.", 
                  "\n")
        }
        if (stand.FRic) {
          hull.all <- convhulln(traits.FRic, "FA")
          FRic.all <- hull.all$vol
        }
      }
    }
  }
  if (!calc.FRic & calc.FDiv) 
    cat("FDiv: Cannot be computed when 'calc.FRic' is FALSE.", 
        "\n")
  if (calc.FRic & calc.FDiv) 
    if (min.nb.sp < 3) 
      if (messages) 
        cat("FDiv: Could not be calculated for communities with <3 functionally singular species.", 
            "\n")
  if (calc.FGR) {
    if (clust.type == "kmeans") {
      tr.clust <- cascadeKM(traits, km.inf.gr, km.sup.gr, 
                            km.iter, km.crit)
      cat("FGR: Summary of kmeans clustering\n")
      cat("\nPartition\n")
      print(tr.clust$partition)
      cat("\nResults\n")
      print(tr.clust$results)
      cat("\nSize\n")
      print(tr.clust$size)
      plot(tr.clust)
      part.names <- colnames(tr.clust$partition)
      part.names <- as.numeric(substr(part.names, 1, 1))
      cat("\nFGR: How many groups?", "\n")
      cut.g <- toupper(scan(file = "", what = "character", 
                            nlines = 1, quiet = T))
      cut.gr <- as.integer(cut.g)
      if (cut.gr < km.inf.gr | cut.gr > km.sup.gr) 
        stop("You must type an integer between 'km.ing.gr' and 'km.sup.gr'.", 
             "\n")
      spfgr.all <- tr.clust$partition[, part.names == cut.gr]
      names(spfgr.all) <- x.rn
    }
    else {
      tr.clust <- hclust(x.dist, method = clust.type)
      plot(tr.clust, main = "Cluster dengrogram of species based on functional traits")
      cat("FGR: Do you want to cut the dendrogram from height or from the number of groups? Type 'h' for height, 'g' for groups.", 
          "\n")
      cut <- toupper(scan(file = "", what = "character", 
                          nlines = 1, quiet = T))
      if (cut == "H") {
        cat("FGR: At what height do you want the dendrogram to be cut?", 
            "\n")
        cut.d <- toupper(scan(file = "", what = "character", 
                              nlines = 1, quiet = T))
        cut.dist <- as.numeric(cut.d)
        spfgr.all <- cutree(tr.clust, h = cut.dist)
      }
      if (cut == "G") {
        cat("FGR: How many groups?", "\n")
        cut.g <- toupper(scan(file = "", what = "character", 
                              nlines = 1, quiet = T))
        cut.gr <- as.integer(cut.g)
        spfgr.all <- cutree(tr.clust, k = cut.gr)
      }
      if (cut != "H" & cut != "G") 
        stop("You must type 'h' or 'g'", "\n")
    }
    a.t <- t(a)
    by.gr <- list(spfgr.all)
    gr.abun <- aggregate(a.t, by.gr, sum)
    lab <- paste("group", gr.abun[, 1], sep = "")
    gr.abun <- data.frame(t(gr.abun[, -1]))
    colnames(gr.abun) <- lab
    rownames(gr.abun) <- rownames(a)
  }
  if (is.matrix(x) | is.data.frame(x) & calc.CWM) {
    CWM <- functcomp(x, a, CWM.type = CWM.type)
  }
  if (calc.CWM & class(x)[1] == "dist" | class(x)[1] == "dissimilarity") 
    if (messages) 
      cat("CWM: When 'x' is a distance matrix, CWM cannot be calculated.", 
          "\n")
  divc <- function(df, dis = NULL, scale = FALSE) {
    if (!inherits(df, "data.frame")) 
      stop("Non convenient df")
    if (any(df < 0)) 
      stop("Negative value in df")
    if (!is.null(dis)) {
      if (!inherits(dis, "dist")) 
        stop("Object of class 'dist' expected for distance")
      dis <- as.matrix(dis)
      if (nrow(df) != nrow(dis)) 
        stop("Non convenient df")
      dis <- as.dist(dis)
    }
    if (is.null(dis)) 
      dis <- as.dist((matrix(1, nrow(df), nrow(df)) - diag(rep(1, 
                                                               nrow(df)))) * sqrt(2))
    div <- as.data.frame(rep(0, ncol(df)))
    names(div) <- "diversity"
    rownames(div) <- names(df)
    for (i in 1:ncol(df)) {
      if (sum(df[, i]) < 1e-16) 
        div[i, ] <- 0
      else div[i, ] <- (t(df[, i]) %*% (as.matrix(dis)^2) %*% 
                          df[, i])/2/(sum(df[, i])^2)
    }
    if (scale == TRUE) {
      divmax <- divcmax(dis)$value
      div <- div/divmax
    }
    return(div)
  }
  #IB added code
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
    AA <- a
    for(i in 1:ncol(a)) AA[,i] <- a[,i]*biomassValue2[i]
    RaoQ <- divc(data.frame(t(AA)), x.dist, scale = scale.RaoQ)    
  }else{    
  RaoQ <- divc(data.frame(t(a)), x.dist, scale = scale.RaoQ) #original line
  }
  ##stop edits
  RaoQ <- RaoQ[, 1]
  names(RaoQ) <- rownames(a)
  disp <- fdisp_w(x.dist, a, Weigthedby = Weigthedby, biomassValue = biomassValue)
  FDis <- disp$FDis
  nbsp <- rep(NA, c)
  names(nbsp) <- row.names(a)
  FRic <- rep(NA, c)
  names(FRic) <- row.names(a)
  FEve <- rep(NA, c)
  names(FEve) <- row.names(a)
  FGR <- rep(NA, c)
  names(FGR) <- row.names(a)
  FDiv <- rep(NA, c)
  names(FDiv) <- row.names(a)
  for (i in 1:c) {
    sppres <- which(a[i, ] > 0)
    S <- length(sppres)
    nbsp[i] <- S
    tr <- data.frame(traits[sppres, ])
    if (calc.FRic) 
      tr.FRic <- data.frame(traits.FRic[sppres, ])
    #start edits IB
    if(Weigthedby != "abundance"){
      ab <- as.matrix(AA[i, sppres])
    }else{
      ab <- as.matrix(a[i, sppres])
    }
    #end
    abundrel <- ab/sum(ab)
    if (calc.FRic) {
      if (all(x.class2 == "factor" | x.class2 == "ordered")) {
        if (length(x.class2) == 1 & x.class2[1] == "ordered") {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
        else {
          if (!stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))
          if (stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))/FRic.all
        }
      }
      else {
        if (dim(tr.FRic)[2] > 1 & nb.sp[i] >= 3) {
          if (warning) 
            thresh <- 4
          if (!warning) 
            thresh <- 3
          if (nb.sp[i] >= thresh) {
            convhull <- convhulln(tr.FRic, "FA")
            if (!stand.FRic) 
              FRic[i] <- convhull$vol
            if (stand.FRic) 
              FRic[i] <- convhull$vol/FRic.all
          }
          else {
          }
        }
        if (dim(tr.FRic)[2] == 1) {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
      }
    }
    if (nb.sp[i] >= 3) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] + 
        abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (m in 1:((S - 1) * S/2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m]/(abund2vect[m])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    }
    if (calc.FDiv & calc.FRic) {
      if (any(x.class2 == "numeric") & dim(tr.FRic)[2] > 
            1 & nb.sp[i] >= 3) {
        vert0 <- convhulln(tr.FRic, "Fx TO 'vert.txt'")
        vert1 <- scan("vert.txt", quiet = T)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        trvertices <- tr.FRic[vertices, ]
        baryv <- apply(trvertices, 2, mean)
        distbaryv <- rep(0, S)
        for (j in 1:S) distbaryv[j] <- (sum((tr.FRic[j, 
                                                     ] - baryv)^2))^0.5
        meandB <- mean(distbaryv)
        devdB <- distbaryv - meandB
        abdev2 <- abundrel * devdB
        ababsdev2 <- abundrel * abs(devdB)
        FDiv[i] <- (sum(abdev2) + meandB)/(sum(ababsdev2) + 
                                             meandB)
      }
    }
    if (calc.FGR) 
      FGR[i] <- length(unique(spfgr.all[sppres]))
  }
  res <- list()
  res$nbsp <- nbsp
  res$sing.sp <- nb.sp
  if (calc.FRic) 
    res$FRic <- FRic
  if (calc.FRic) 
    res$qual.FRic <- qual.FRic
  res$FEve <- FEve
  if (calc.FDiv) 
    res$FDiv <- FDiv
  res$FDis <- FDis
  res$RaoQ <- RaoQ
  if (calc.FGR) {
    res$FGR <- FGR
    res$spfgr <- spfgr.all
    res$gr.abun <- gr.abun
  }
  if (is.matrix(x) | is.data.frame(x) & calc.CWM) 
    res$CWM <- CWM
  if (print.pco) {
    res$x.values <- x.pco$eig
    res$x.axes <- x.pco$li
  }
  invisible(res)
}