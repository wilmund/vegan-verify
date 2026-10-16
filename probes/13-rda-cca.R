# rda() / cca(): eigenvalues and inertia recomputed from the definitional
# linear algebra in base R.
#
# rda: the init transform is Y = center(comm) / sqrt(n-1) (optionally
# scale()d first), so unconstrained eigenvalues are svd(Y)$d^2 and total
# inertia is sum(Y^2) (= sum of column variances). cca: the init transform
# is the chi-square standardized Q = (P - rc') / sqrt(rc') with P = comm /
# grand total, r/c the margin sums of P; eigenvalues svd(Q)$d^2, inertia
# sum(Q^2). Constrained case: project onto the (row-weighted, for cca)
# centered constraint matrix by QR least squares — constrained eigenvalues
# are svd of the fitted part, residual ("unconstrained") eigenvalues svd of
# the residual, and total = constrained + unconstrained inertia must hold
# as an accounting identity. inertcomp() splits inertia by species, so its
# grand sum must return the total. Zero eigenvalues are dropped with
# vegan's own threshold (max(ZERO, ZERO * lambda[1]), ZERO = sqrt(eps)) so
# both sides report the same length.

comm <- as.matrix(read.csv("data/community6.csv")) * 1.0
env <- as.data.frame(read.csv("data/env6.csv"))
n <- nrow(comm)
ZERO <- sqrt(.Machine$double.eps)
dropz <- function(l) l[abs(l) > max(ZERO, ZERO * l[1])]

## ---- rda, unconstrained -------------------------------------------------
Yr <- scale(comm, scale = FALSE) / sqrt(n - 1)
f <- rda(comm)
check_close("rda unconstrained eigenvalues = svd(centered Y / sqrt(n-1))^2",
            dropz(svd(Yr)$d^2), unname(f$CA$eig))
check_close("rda total inertia = sum(Y^2) = sum of column variances",
            sum(Yr^2), f$tot.chi)
check_close("rda total inertia (cross-check vs apply(var))",
            sum(apply(comm, 2, var)), f$tot.chi)

## rda with scale = TRUE: same construction on the scaled matrix
Ys <- scale(comm) / sqrt(n - 1)
fs <- rda(comm, scale = TRUE)
check_close("rda scale=TRUE eigenvalues = svd(scaled Y / sqrt(n-1))^2",
            dropz(svd(Ys)$d^2), unname(fs$CA$eig))

## ---- rda, constrained ---------------------------------------------------
X <- scale(as.matrix(env[, c("e1", "e2")]), scale = FALSE)
QR <- qr(X)
Yfit <- qr.fitted(QR, Yr)
f2 <- rda(comm ~ e1 + e2, data = env)
check_close("rda constrained eigenvalues = svd(fitted part)^2",
            dropz(svd(Yfit)$d^2), unname(f2$CCA$eig))
check_close("rda residual eigenvalues = svd(Y - fitted)^2",
            dropz(svd(Yr - Yfit)$d^2), unname(f2$CA$eig))
check_close("rda constrained inertia = sum(fitted^2)",
            sum(Yfit^2), f2$CCA$tot.chi)
check_close("rda inertia accounting: total = constrained + unconstrained",
            f2$tot.chi, f2$CCA$tot.chi + f2$CA$tot.chi)
check_close("rda inertcomp grand sum = total inertia",
            sum(inertcomp(f2)), f2$tot.chi)

## ---- cca, unconstrained -------------------------------------------------
P <- comm / sum(comm)
rw <- rowSums(P)
cw <- colSums(P)
Qb <- (P - outer(rw, cw)) / sqrt(outer(rw, cw))
g <- cca(comm)
check_close("cca unconstrained eigenvalues = svd(chi-standardized Q)^2",
            dropz(svd(Qb)$d^2), unname(g$CA$eig))
check_close("cca total inertia = sum(Q^2)", sum(Qb^2), g$tot.chi)

## ---- cca, constrained (row-weighted projection) -------------------------
Xc <- as.matrix(env[, c("e1", "e2")])
Xc <- scale(Xc, center = apply(Xc, 2, weighted.mean, w = rw), scale = FALSE)
Xc <- Xc * sqrt(rw)
Qfit <- qr.fitted(qr(Xc), Qb)
g2 <- cca(comm ~ e1 + e2, data = env)
check_close("cca constrained eigenvalues = svd(weighted fitted part)^2",
            dropz(svd(Qfit)$d^2), unname(g2$CCA$eig))
check_close("cca residual eigenvalues = svd(Q - fitted)^2",
            dropz(svd(Qb - Qfit)$d^2), unname(g2$CA$eig))
check_close("cca constrained inertia = sum(fitted^2)",
            sum(Qfit^2), g2$CCA$tot.chi)
check_close("cca inertia accounting: total = constrained + unconstrained",
            g2$tot.chi, g2$CCA$tot.chi + g2$CA$tot.chi)
check_close("cca inertcomp grand sum = total inertia",
            sum(inertcomp(g2)), g2$tot.chi)
