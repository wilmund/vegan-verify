#!/usr/bin/env Rscript
# run-probes.R — run every probe in probes/ against the installed vegan.
#
# Usage: Rscript run-probes.R [label]
# Prints one PASS/FAIL line per check and a summary; exit status 1 if any
# check failed. The output is the body of a dated report in reports/.

label <- commandArgs(trailingOnly = TRUE)
label <- if (length(label)) label[[1]] else "unlabelled"

source("harness.R")
suppressPackageStartupMessages(library(vegan))

cat(sprintf("vegan-verify run: %s\n", label))
cat(sprintf("R %s | vegan %s | loaded from %s | %s\n\n", getRversion(),
            as.character(packageVersion("vegan")),
            dirname(system.file(package = "vegan")),
            format(Sys.time(), "%Y-%m-%d %H:%M %Z")))

for (f in sort(list.files("probes", pattern = "\\.R$", full.names = TRUE))) {
  cat("--", basename(f), "--\n")
  sys.source(f, envir = new.env(parent = globalenv()))
  cat("\n")
}

fails <- harness_summary()
quit(save = "no", status = if (fails > 0) 1L else 0L)
