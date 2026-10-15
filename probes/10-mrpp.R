# mrpp(): delta, E(delta), A and the exact permutation p recomputed from
# their definitions on a hand-built Bray-Curtis matrix. delta = weighted
# mean of within-group mean dissimilarities; the documented weight types
# are w = n_g (1), n_g - 1 (2), n_g(n_g-1)/2 (3) -- the grouping here is
# UNBALANCED (2 + 4) precisely so the three weightings give different
# numbers. Under complete enumeration every pair is equally likely to be
# within-group, so E(delta) = the plain mean of all pairwise
# dissimilarities (closed form, no permutations needed) and A = 1 -
# delta/E(delta). Exact p = share of the 720 orderings with delta* <=
# delta_obs (small delta = concentration); vegan gets the 719 non-identity
# orderings so Pvalue = (hits+1)/720 must equal the exact p.

comm <- as.matrix(read.csv("data/community6.csv"))
grp <- factor(c("A", "A", "B", "B", "B", "B"))
n <- nrow(comm)
EPS <- sqrt(.Machine$double.eps)

bray <- function(x, y) sum(abs(x - y)) / sum(x + y)
D <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) D[i, j] <- bray(comm[i, ], comm[j, ])

classmean <- function(g) {
  vapply(levels(g), function(lev) {
    idx <- which(g == lev)
    pairs <- t(combn(idx, 2))
    mean(D[pairs])
  }, numeric(1))
}
sizes <- as.vector(table(grp))
wts <- list(sizes, sizes - 1, sizes * (sizes - 1) / 2)
delta_of <- function(g, w) {
  cm <- classmean(g)
  sum(w * cm) / sum(w)
}

cm_obs <- classmean(grp)
delta_obs <- delta_of(grp, wts[[1]])
E_delta <- mean(D[lower.tri(D)])
A_obs <- 1 - delta_obs / E_delta

P <- all_perms(n)
ddist <- apply(P, 1L, function(p) delta_of(grp[p], wts[[1]]))
p_exact <- sum(ddist <= delta_obs + EPS) / nrow(P)

pm <- perms_no_identity(n)
fit <- mrpp(as.dist(D), grp, permutations = pm, weight.type = 1)

check_close("mrpp delta (weight.type 1)", delta_obs, fit$delta)
check_close("mrpp classdelta per group", as.vector(cm_obs),
            as.vector(fit$classdelta))
check_close("mrpp E(delta) = mean of all pairwise dissimilarities",
            E_delta, fit$E.delta)
check_close("mrpp A = 1 - delta/E(delta)", A_obs, fit$A)
check_close("mrpp exact p over all 720 orderings", p_exact, fit$Pvalue)

# Weight types 2 and 3 reweight delta; the unbalanced grouping makes them
# distinct numbers, each with its own exact enumeration p.
for (wt in 2:3) {
  dobs <- delta_of(grp, wts[[wt]])
  dd <- apply(P, 1L, function(p) delta_of(grp[p], wts[[wt]]))
  fitw <- mrpp(as.dist(D), grp, permutations = pm, weight.type = wt)
  check_close(sprintf("mrpp delta (weight.type %d)", wt), dobs, fitw$delta)
  check_close(sprintf("mrpp exact p (weight.type %d)", wt),
              sum(dd <= dobs + EPS) / nrow(P), fitw$Pvalue)
}
