#!/bin/sh
# run-suite.sh — run the probe suite against BOTH the CRAN release and the
# vegan dev head, writing dated reports to reports/. Expects to run inside
# a rocker/r2u container with CRAN vegan installed system-wide and a fresh
# vegan clone at /build (cp -r the read-only mount first; installing from
# a ro bind mount fails).
#
#   docker run -d --name vegan-probe -v /tmp/vegan:/repo:ro rocker/r2u:latest sleep infinity
#   docker exec vegan-probe Rscript -e 'install.packages("vegan")'
#   docker cp . vegan-probe:/vv
#   docker exec vegan-probe sh -c 'cp -r /repo /build && mkdir -p /tmp/devlib && R CMD INSTALL -l /tmp/devlib /build'
#   docker exec -w /vv vegan-probe sh ./run-suite.sh <dev-sha>
#
# Exit status: number of failing targets (0 = both green).

set -u
DEVSHA="${1:-unknown-sha}"
# Reports are append-only records: a UTC date+time stamp in the filename so
# a same-day re-run can never silently overwrite an earlier report.
STAMP=$(date -u +%Y-%m-%d-%H%MZ)
mkdir -p reports
FAILS=0

Rscript run-probes.R "CRAN release" > "reports/$STAMP-cran.txt" 2>&1 || FAILS=$((FAILS+1))
R_LIBS=/tmp/devlib Rscript run-probes.R "dev head $DEVSHA" > "reports/$STAMP-dev-$DEVSHA.txt" 2>&1 || FAILS=$((FAILS+1))

tail -1 "reports/$STAMP-cran.txt"
tail -1 "reports/$STAMP-dev-$DEVSHA.txt"
exit $FAILS
