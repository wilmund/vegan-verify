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
