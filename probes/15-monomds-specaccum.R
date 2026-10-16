# monoMDS(): not the configuration (stochastic optimization), but the
# STRESS OF THE RETURNED CONFIGURATION, recomputed from Kruskal's
# definition: dhat = monotone (isotonic) regression of the configuration
# distances on the dissimilarity order — implemented here as a plain
# pool-adjacent-violators loop in base R — and
#   stress 1 = sqrt(sum((d - dhat)^2) / sum(d^2))
#   stress 2 = sqrt(sum((d - dhat)^2) / sum((d - mean(d))^2))
# model = "linear" replaces the monotone fit with ordinary least squares.
# vegan's Fortran returns $dist and $dhat sorted by ascending
# dissimilarity on the MONOTONE path (MONREG requires sorted input), but in
# the original lower-triangle order on the linear path (LINREG does not
# sort) — the correspondence checks use the matching order for each. The
# dissimilarity matrix is a fixed hand-written one whose rank structure
# deliberately conflicts with any 1-D embedding: at k = 1 the fit is bad
# (stress ~0.4), the isotonic regression must pool violators, and the
# stress formula is exercised far from the trivial dhat = d case. The
# start configuration is passed explicitly, so no RNG is involved. All 15
# dissimilarities are distinct (asserted), so tie handling is moot.
#
# specaccum(method = "exact"): the expected species accumulation curve has
# a closed form, E[S_k] = sum_i (1 - choose(n - f_i, k) / choose(n, k))
# with f_i = number of sites holding species i — and both it and
# specaccum(method = "random") handed ALL n! site orderings must equal the
# brute-force average over every ordering. The random-method sd must equal
# the brute-force sd over all orderings. (The exact-method sd is Kindt's
# conditioned approximation, not the permutation sd — deliberately not
# compared against this oracle.)

env <- as.matrix(read.csv("data/env6.csv"))

## ---- monoMDS ------------------------------------------------------------
vals <- c(0.71, 0.23, 0.58, 0.91, 0.34,
                0.47, 0.12, 0.66, 0.82,
                      0.95, 0.28, 0.53,
                            0.39, 0.77,
                                  0.19)
D <- matrix(0, 6, 6)
D[lower.tri(D)] <- vals
D <- D + t(D)
diss <- D[lower.tri(D)]
check_true("monoMDS input dissimilarities are tie-free",
           length(unique(diss)) == length(diss))

## pool-adjacent-violators: isotonic regression with unit weights,
## y already ordered by ascending dissimilarity
pava <- function(y) {
  nb <- length(y); val <- y; w <- rep(1, nb); k <- nb
  repeat {
    viol <- which(diff(val[1:k]) < 0)
    if (!length(viol)) break
    i <- viol[1]
    val[i] <- (w[i] * val[i] + w[i + 1] * val[i + 1]) / (w[i] + w[i + 1])
    w[i] <- w[i] + w[i + 1]
    if (i + 1 < k) {
      val[(i + 1):(k - 1)] <- val[(i + 2):k]
      w[(i + 1):(k - 1)] <- w[(i + 2):k]
    }
    k <- k - 1
  }
  rep(val[1:k], times = w[1:k])
}

ord <- order(diss)
y0 <- scale(env[, 1, drop = FALSE], scale = FALSE)  # deterministic start

for (ss in 1:2) {
  m <- monoMDS(as.dist(D), y = y0, k = 1, model = "global", stress = ss)
  dc <- as.matrix(dist(m$points))
  dconf <- dc[lower.tri(dc)]
  dhat <- pava(dconf[ord])
  s <- if (ss == 1) sqrt(sum((dconf[ord] - dhat)^2) / sum(dconf^2))
       else sqrt(sum((dconf[ord] - dhat)^2) / sum((dconf - mean(dconf))^2))
  check_true(sprintf("monoMDS stress-%d run is a non-trivial fit (stress > 0.1)", ss),
             m$stress > 0.1)
  check_close(sprintf("monoMDS $dist (form %d) = distances of returned $points", ss),
              dconf[ord], m$dist)
  check_close(sprintf("monoMDS $dhat (form %d) = base-R isotonic regression", ss),
              dhat, m$dhat)
  check_close(sprintf("monoMDS stress-%d = Kruskal formula on returned configuration", ss),
              s, m$stress)
}

## model = "linear": dhat is an ordinary least-squares fit
ml <- monoMDS(as.dist(D), y = y0, k = 1, model = "linear")
dcl <- as.matrix(dist(ml$points))
dconfl <- dcl[lower.tri(dcl)]
lmfit <- lm(dconfl ~ diss)
check_close("monoMDS linear $dist = distances of returned $points (original order)",
            dconfl, ml$dist)
check_close("monoMDS linear $dhat = lm fitted values (original order)",
            unname(fitted(lmfit)), ml$dhat)
check_close("monoMDS linear stress = formula with lm residuals",
            sqrt(sum(residuals(lmfit)^2) / sum(dconfl^2)), ml$stress)

## ---- specaccum ----------------------------------------------------------
comm5 <- as.matrix(read.csv("data/community.csv"))
n5 <- nrow(comm5)
P5 <- all_perms(n5)
acc <- apply(P5, 1, function(p)
  rowSums(apply(comm5[p, , drop = FALSE], 2, cumsum) > 0))
bf_mean <- rowMeans(acc)
bf_sd <- apply(acc, 1, sd)

sx <- specaccum(comm5, method = "exact")
check_close("specaccum exact = brute-force mean over all 120 orderings",
            bf_mean, sx$richness)
freq <- colSums(comm5 > 0)
check_close("specaccum exact = closed form sum(1 - choose(n-f,k)/choose(n,k))",
            sapply(1:n5, function(k) sum(1 - choose(n5 - freq, k) / choose(n5, k))),
            unname(sx$richness))

sr <- specaccum(comm5, method = "random", permutations = P5)
check_close("specaccum random over all orderings: mean = brute force",
            bf_mean, sr$richness)
check_close("specaccum random over all orderings: sd = brute force",
            bf_sd, sr$sd)
