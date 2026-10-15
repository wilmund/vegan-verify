# fisher.alpha(): verified by solution-substitution — the returned alpha
# must satisfy S = alpha * log(1 + N/alpha). No reimplementation needed.
# Tolerance 1e-5 (relative): alpha is a numerical root (nlm inside vegan),
# and observed substitution residuals on these data are up to ~2e-6 — the
# solver's convergence error, not a defect. A wrong equation or mishandled
# input produces residuals orders of magnitude larger.

comm <- as.matrix(read.csv("data/community.csv"))
for (i in seq_len(nrow(comm))) {
  x <- comm[i, ]
  S <- sum(x > 0); N <- sum(x)
  a <- as.numeric(fisher.alpha(x))
  check_close(sprintf("fisher.alpha substitution row %d", i),
              S, a * log(1 + N / a), tol = 1e-5)
}
