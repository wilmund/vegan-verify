# anosim(): the R statistic recomputed from its definition -- rank all
# n(n-1)/2 dissimilarities (average ranks on ties, the standard
# convention), then R = (mean rank between-group - mean rank within-group)
# divided by M/2 with M = n(n-1)/2 -- on a hand-built Bray-Curtis matrix.
# Exact p-value by enumerating all 720 label orderings; vegan gets the 719
# non-identity orderings so signif = (hits+1)/720 must equal the exact p.

comm <- as.matrix(read.csv("data/community6.csv"))
grp <- factor(c("A", "A", "A", "B", "B", "B"))
n <- nrow(comm)
EPS <- sqrt(.Machine$double.eps)

bray <- function(x, y) sum(abs(x - y)) / sum(x + y)
D <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) D[i, j] <- bray(comm[i, ], comm[j, ])

lt <- lower.tri(D)
rk <- rank(D[lt])          # base-R rank, average ties
M <- sum(lt)
div <- M / 2

R_of <- function(g) {
  within <- outer(g, g, "==")[lt]
  mean(rk[!within]) / div - mean(rk[within]) / div
}

R_obs <- R_of(grp)
P <- all_perms(n)
Rdist <- apply(P, 1L, function(p) R_of(grp[p]))
p_exact <- sum(Rdist >= R_obs - EPS) / nrow(P)

fit <- anosim(as.dist(D), grp, permutations = perms_no_identity(n))

check_close("anosim R statistic from rank definition", R_obs, fit$statistic)
check_close("anosim exact p over all 720 orderings", p_exact, fit$signif)
