# diversity(): Shannon, Simpson and inverse Simpson recomputed from their
# definitions with a plain loop over proportions. Edge cases with known
# closed-form answers: a single-species community and an all-equal one.

comm <- as.matrix(read.csv("data/community.csv"))

shannon_hand    <- function(x) { p <- x[x > 0] / sum(x); -sum(p * log(p)) }
simpson_hand    <- function(x) { p <- x / sum(x); 1 - sum(p^2) }
invsimpson_hand <- function(x) { p <- x / sum(x); 1 / sum(p^2) }

for (i in seq_len(nrow(comm))) {
  x <- comm[i, ]
  check_close(sprintf("diversity shannon row %d", i),
              shannon_hand(x), unname(diversity(x, "shannon")))
  check_close(sprintf("diversity simpson row %d", i),
              simpson_hand(x), unname(diversity(x, "simpson")))
  check_close(sprintf("diversity invsimpson row %d", i),
              invsimpson_hand(x), unname(diversity(x, "invsimpson")))
}

one_sp <- c(7, 0, 0)
check_close("diversity shannon single-species = 0", 0, unname(diversity(one_sp, "shannon")))
check_close("diversity simpson single-species = 0", 0, unname(diversity(one_sp, "simpson")))
check_close("diversity invsimpson single-species = 1", 1, unname(diversity(one_sp, "invsimpson")))

eq <- rep(4, 6)
check_close("diversity shannon all-equal = log(S)", log(6), unname(diversity(eq, "shannon")))
check_close("diversity invsimpson all-equal = S", 6, unname(diversity(eq, "invsimpson")))
