# mantel(): the statistic is the plain product-moment correlation of the
# two lower triangles, recomputed here from the explicit Pearson formula
# (and, for method = "spearman", the same formula on base-R ranks) on
# hand-built matrices: Bray-Curtis on the community, Euclidean on the
# environment. Exact p by enumerating all 720 simultaneous row/column
# orderings of one matrix; over the COMPLETE group it does not matter
# which of the two matrices is permuted (the two distributions are images
# of each other under p -> p^-1). vegan gets the 719 non-identity
# orderings so signif = (hits+1)/720 must equal the exact p.

comm <- as.matrix(read.csv("data/community6.csv"))
env <- as.matrix(read.csv("data/env6.csv"))
n <- nrow(comm)
EPS <- sqrt(.Machine$double.eps)

bray <- function(x, y) sum(abs(x - y)) / sum(x + y)
euc <- function(x, y) sqrt(sum((x - y)^2))
D1 <- matrix(0, n, n)
D2 <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) {
  D1[i, j] <- bray(comm[i, ], comm[j, ])
  D2[i, j] <- euc(env[i, ], env[j, ])
}

lt <- lower.tri(D1)
pearson <- function(x, y) {
  xm <- x - mean(x); ym <- y - mean(y)
  sum(xm * ym) / sqrt(sum(xm^2) * sum(ym^2))
}

x <- D1[lt]
P <- all_perms(n)
pm <- perms_no_identity(n)

for (meth in c("pearson", "spearman")) {
  f <- if (meth == "pearson") identity else rank
  stat_of <- function(p) {
    D2p <- D2[p, p]
    pearson(f(x), f(D2p[lt]))
  }
  r_obs <- stat_of(seq_len(n))
  rdist <- apply(P, 1L, stat_of)
  p_exact <- sum(rdist >= r_obs - EPS) / nrow(P)

  fit <- mantel(as.dist(D1), as.dist(D2), method = meth, permutations = pm)
  check_close(sprintf("mantel %s statistic = cor of lower triangles", meth),
              r_obs, fit$statistic)
  check_close(sprintf("mantel %s exact p over all 720 orderings", meth),
              p_exact, fit$signif)
}
