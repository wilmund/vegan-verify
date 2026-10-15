# harness.R — minimal check framework sourced by run-probes.R.
#
# Probes call check_identical() for integer/structural results and
# check_close() for floating-point results. Default tolerance is a
# RELATIVE 1e-12: each oracle recomputes its quantity in a different
# operation order than vegan, so bit-identity is not guaranteed, but any
# real defect is many orders of magnitude above 1e-12. Probes that need
# a wider tolerance must justify it in a comment at the call site.

RESULTS <- new.env(parent = emptyenv())
RESULTS$rows <- list()

.record <- function(name, ok, detail = "") {
  RESULTS$rows[[length(RESULTS$rows) + 1L]] <- list(name = name, ok = ok, detail = detail)
  cat(sprintf("%s  %s%s\n", if (ok) "PASS" else "FAIL", name,
              if (ok || !nzchar(detail)) "" else paste0("  [", detail, "]")))
  invisible(ok)
}

check_identical <- function(name, expected, got) {
  ok <- identical(expected, got)
  .record(name, ok, if (!ok) paste0(
    "expected ", paste(format(expected), collapse = " "),
    " | got ", paste(format(got), collapse = " ")) else "")
}

check_close <- function(name, expected, got, tol = 1e-12) {
  expected <- as.vector(expected); got <- as.vector(got)
  cmp <- all.equal(expected, got, tolerance = tol)
  ok <- isTRUE(cmp) && length(expected) == length(got)
  .record(name, ok, if (!ok) paste0(
    "expected ", paste(format(expected, digits = 17), collapse = " "),
    " | got ", paste(format(got, digits = 17), collapse = " "),
    " | ", paste(as.character(cmp), collapse = "; ")) else "")
}

check_true <- function(name, cond, detail = "") {
  .record(name, isTRUE(cond), detail)
}

harness_summary <- function() {
  ok <- vapply(RESULTS$rows, function(r) r$ok, logical(1))
  cat(sprintf("\nSUMMARY: %d checks, %d passed, %d FAILED\n",
              length(ok), sum(ok), sum(!ok)))
  if (any(!ok)) {
    cat("Failed checks:\n")
    for (r in RESULTS$rows) if (!r$ok) cat("  -", r$name, "\n")
  }
  sum(!ok)
}

# --- Tier 2 helpers: exhaustive permutation enumeration ------------------
#
# all_perms(n) returns every ordering of 1..n as an n! x n integer matrix,
# built by plain recursion in base R (no permute package: the enumeration
# is part of the oracle, so it must not come from vegan's dependency).
# perms_no_identity(n) drops the identity row: vegan's permutation tests
# compute p = (hits + 1)/(nperm + 1) over a SUPPLIED matrix, so feeding the
# n! - 1 non-identity orderings makes the reported p exactly equal to the
# full-enumeration p (the +1 is the identity, which always "hits").

all_perms <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  sub <- all_perms(n - 1L)
  rows <- nrow(sub)
  out <- matrix(0L, rows * n, n)
  for (k in seq_len(n)) {
    rest <- setdiff(seq_len(n), k)
    out[(k - 1L) * rows + seq_len(rows), ] <-
      cbind(rep(k, rows), matrix(rest[sub], rows, n - 1L))
  }
  out
}

perms_no_identity <- function(n) {
  P <- all_perms(n)
  P[!apply(P, 1L, function(r) all(r == seq_len(n))), , drop = FALSE]
}
