# rarefy(): the expected species count E[S_n] computed two independent
# ways — the closed hypergeometric formula, and brute-force enumeration of
# EVERY size-n subsample of a 10-individual toy community (the strongest
# oracle available: the answer is the average over all possible draws).

x <- c(spA = 4, spB = 3, spC = 2, spD = 1)   # N = 10 individuals
N <- sum(x)
inds <- rep(names(x), x)

for (n in c(2, 5, 8)) {
  hyper <- sum(1 - choose(N - x, n) / choose(N, n))
  subs <- combn(N, n)
  brute <- mean(apply(subs, 2, function(idx) length(unique(inds[idx]))))
  check_close(sprintf("rarefy hypergeometric = brute force, n=%d", n), brute, hyper)
  check_close(sprintf("rarefy vegan = formula, n=%d", n), hyper, as.vector(rarefy(x, n)))
}

# Whole community matrix at its minimum row total.
comm <- as.matrix(read.csv("data/community.csv"))
nmin <- min(rowSums(comm))
hyper_rows <- apply(comm, 1, function(r) sum(1 - choose(sum(r) - r, nmin) / choose(sum(r), nmin)))
check_close("rarefy matrix rows vs formula", unname(hyper_rows), as.vector(rarefy(comm, nmin)))
