# betadisper(): distances to group centroids recomputed from first
# principles (Anderson 2006 / stats::cmdscale PCoA construction), including
# the negative-eigenvalue correction — the subtle spot where a regression
# would hide. Bray-Curtis on count data is non-Euclidean, so the PCoA here
# genuinely has a negative eigenvalue (asserted below, otherwise the
# correction path would be untested dead weight).
#
# Hand construction: square the dissimilarities, double-center (-x/2),
# eigen-decompose; drop near-zero eigenvalues with vegan's own threshold;
# scale eigenvectors by sqrt(|eigenvalue|). Squared distance to centroid =
# (squared distance on real/positive axes) MINUS (squared distance on
# imaginary/negative axes); zij = sqrt of that. type = "centroid" so the
# centroid is a plain group mean (type = "median" is an LAD optimization,
# an algorithm rather than a definition — out of scope for this oracle).
# bias.adjust multiplies zij by sqrt(n_g / (n_g - 1)).
#
# The degenerate corner where dist.neg > dist.pos (complex distance clamped
# to zero, vegan GitHub issue #306) is not reachable with well-separated
# groups like these; it is not exercised here.

comm <- as.matrix(read.csv("data/community6.csv")) * 1.0
n <- nrow(comm)
ZERO <- sqrt(.Machine$double.eps)

## hand-built Bray-Curtis distances (double loop, no vegdist)
D <- matrix(0, n, n)
for (i in 1:n) for (j in 1:n)
  D[i, j] <- sum(abs(comm[i, ] - comm[j, ])) / sum(comm[i, ] + comm[j, ])
d <- as.dist(D)

## PCoA by hand
x <- matrix(0, n, n)
x[row(x) > col(x)] <- as.vector(d)^2
x <- x + t(x)
x <- sweep(x, 2, colMeans(x))
x <- sweep(x, 1, rowMeans(x))
e <- eigen(-x / 2, symmetric = TRUE)
keep <- abs(e$values) > max(ZERO, ZERO * e$values[1])
eig <- e$values[keep]
V <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(abs(eig)), length(eig))
pos <- eig > 0

check_true("PCoA of Bray on counts has a negative eigenvalue (correction path exercised)",
           any(!pos))

## distances to group centroids with the negative-axis correction
grp <- factor(c("a", "a", "a", "b", "b", "b"))
cent <- apply(V, 2, function(col) tapply(col, grp, mean))
dp <- rowSums((V[, pos, drop = FALSE] - cent[grp, pos, drop = FALSE])^2)
dn <- rowSums((V[, !pos, drop = FALSE] - cent[grp, !pos, drop = FALSE])^2)
zij <- sqrt(dp - dn)

bd <- betadisper(d, grp, type = "centroid")
check_close("betadisper eigenvalues = hand PCoA eigenvalues",
            eig, unname(bd$eig))
check_close("betadisper distances = sqrt(pos-axis sq dist - neg-axis sq dist)",
            zij, unname(bd$distances))
check_close("betadisper group.distances = group means of zij",
            as.vector(tapply(zij, grp, mean)), as.vector(bd$group.distances))

bd2 <- betadisper(d, grp, type = "centroid", bias.adjust = TRUE)
ng <- as.vector(table(grp))
check_close("betadisper bias.adjust multiplies by sqrt(n/(n-1))",
            zij * sqrt(ng[grp] / (ng[grp] - 1)), unname(bd2$distances))
