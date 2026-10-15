# vegdist(): Bray-Curtis, quantitative Jaccard, Euclidean, Kulczynski and
# Gower recomputed pairwise from their textbook formulas. Jaccard is
# computed via the independent Ruzicka route 1 - sum(min)/sum(max), which
# equals 2B/(1+B) algebraically but shares no code path with vegan's.
# Gower here = mean over columns of |x-y| / column range; requires every
# column range to be non-degenerate (asserted below).

comm <- as.matrix(read.csv("data/community.csv"))
pairs <- t(combn(nrow(comm), 2))
rng <- apply(comm, 2, function(col) diff(range(col)))
stopifnot(all(rng > 0))

hand <- function(method, x, y) {
  switch(method,
    bray       = sum(abs(x - y)) / sum(x + y),
    jaccard    = 1 - sum(pmin(x, y)) / sum(pmax(x, y)),
    euclidean  = sqrt(sum((x - y)^2)),
    kulczynski = 1 - 0.5 * (sum(pmin(x, y)) / sum(x) + sum(pmin(x, y)) / sum(y)),
    gower      = mean(abs(x - y) / rng))
}

for (m in c("bray", "jaccard", "euclidean", "kulczynski", "gower")) {
  d <- as.matrix(vegdist(comm, method = m))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    check_close(sprintf("vegdist %s (%d,%d)", m, i, j),
                hand(m, comm[i, ], comm[j, ]), d[i, j])
  }
}

# Internal-consistency triangulation: the designdist() route to Bray-Curtis
# must agree with vegdist(). (Both are vegan, but independent code paths.)
dd <- as.matrix(designdist(comm, "(A+B-2*J)/(A+B)", terms = "minimum"))
db <- as.matrix(vegdist(comm, "bray"))
check_close("designdist minimum-terms route = vegdist bray",
            db[lower.tri(db)], dd[lower.tri(dd)])
