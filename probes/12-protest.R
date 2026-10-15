# protest(): the Procrustes correlation recomputed from the closed-form
# solution of the orthogonal Procrustes problem (Gower & Dijksterhuis):
# after centering both configurations and scaling each to unit sum of
# squares, the minimized symmetric m^2 is 1 - (sum of singular values of
# X'Y)^2, and protest's statistic t0 = sqrt(1 - m^2) = that sum of
# singular values. svd() here is base R, not vegan. Exact p by enumerating
# all 720 row orderings of Y; vegan gets the 719 non-identity orderings so
# signif = (hits+1)/720 must equal the exact p.

comm <- as.matrix(read.csv("data/community6.csv"))
env <- as.matrix(read.csv("data/env6.csv"))
n <- nrow(comm)
EPS <- sqrt(.Machine$double.eps)

X <- env[, 1:2]
Y <- comm[, c(1, 3)] * 1.0

unitize <- function(M) {
  M <- scale(M, scale = FALSE)
  M / sqrt(sum(M^2))
}
Xc <- unitize(X)
Yc <- unitize(Y)

t_of <- function(p) sum(svd(crossprod(Xc, Yc[p, ]))$d)

t0 <- t_of(seq_len(n))
P <- all_perms(n)
tdist <- apply(P, 1L, t_of)
p_exact <- sum(tdist >= t0 - EPS) / nrow(P)

fit <- protest(X, Y, permutations = perms_no_identity(n))

check_close("protest t0 = sum of singular values of X'Y", t0, fit$t0)
check_close("protest ss = 1 - t0^2", 1 - t0^2, fit$ss)
check_close("protest exact p over all 720 orderings", p_exact, fit$signif)
