# adonis2(): one-factor PERMANOVA verified by exhaustive enumeration.
# The oracle recomputes the pseudo-F from the McArdle-Anderson sums of
# squares on a hand-built Bray-Curtis matrix (base-R double loop, no
# vegdist): SS_total = sum over all pairs of d^2 / n, SS_within = sum over
# groups of within-group pair d^2 / n_g, F = (SS_b/(a-1)) / (SS_w/(n-a)).
# It then enumerates all 6! = 720 site orderings and computes the exact
# permutation p-value. vegan is handed the 719 non-identity orderings
# explicitly, so its Pr(>F) = (hits+1)/720 must equal the exact p. The tie
# tolerance EPS = sqrt(.Machine$double.eps) mirrors the convention vegan's
# permutation tests document for counting equal statistics as hits.

comm <- as.matrix(read.csv("data/community6.csv"))
grp <- factor(c("A", "A", "A", "B", "B", "B"))
n <- nrow(comm)
a <- nlevels(grp)
EPS <- sqrt(.Machine$double.eps)

bray <- function(x, y) sum(abs(x - y)) / sum(x + y)
D <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) D[i, j] <- bray(comm[i, ], comm[j, ])

pseudoF <- function(g) {
  SSt <- 0
  for (i in 1:(n - 1)) for (j in (i + 1):n) SSt <- SSt + D[i, j]^2
  SSt <- SSt / n
  SSw <- 0
  for (lev in levels(g)) {
    idx <- which(g == lev)
    s <- 0
    for (ii in seq_along(idx)) for (jj in seq_along(idx))
      if (ii < jj) s <- s + D[idx[ii], idx[jj]]^2
    SSw <- SSw + s / length(idx)
  }
  SSb <- SSt - SSw
  list(F = (SSb / (a - 1)) / (SSw / (n - a)),
       R2 = SSb / SSt, SSb = SSb, SSw = SSw, SSt = SSt)
}

obs <- pseudoF(grp)
P <- all_perms(n)
Fdist <- apply(P, 1L, function(p) pseudoF(grp[p])$F)
p_exact <- sum(Fdist >= obs$F - EPS) / nrow(P)

pm <- perms_no_identity(n)
fit <- adonis2(as.dist(D) ~ grp, data = data.frame(grp = grp),
               permutations = pm)

check_close("adonis2 pseudo-F = hand McArdle-Anderson SS", obs$F, fit$F[1])
check_close("adonis2 R2 = SS_between/SS_total", obs$R2, fit$R2[1])
check_close("adonis2 SumOfSqs (model, residual, total)",
            c(obs$SSb, obs$SSw, obs$SSt), fit$SumOfSqs)
check_close("adonis2 exact p over all 720 orderings",
            p_exact, fit$`Pr(>F)`[1])
