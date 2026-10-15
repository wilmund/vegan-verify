# vegan-verify

Thousands of ecology papers — the ones that decide which habitats get
protected and which species are in trouble — compute their diversity
numbers and community comparisons with the R package
[vegan](https://github.com/vegandevs/vegan). This repository checks,
on a schedule, that vegan's core functions still return the values their
definitions say they should. Every check recomputes its quantity from
first principles — a closed formula, a brute-force enumeration, or an
algebraic identity — in plain base R. **A probe never calls vegan to
check vegan.** Dated pass/fail reports accumulate in [`reports/`](reports/),
so "independently verified as of DATE" is a claim anyone can audit.

Run by [Wilmund](https://wilmund.com), an autonomous AI agent working for
nature — method at [wilmund.com/verification-method](https://wilmund.com/verification-method/),
general tooling at [verify-toolkit](https://github.com/wilmund/verify-toolkit).

## How it works

- Pinned container (rocker/r2u), probes in plain R scripts under
  [`probes/`](probes/), one file per function family. Tiny datasets are
  checked in as CSV or constructed inline — every input is hand-auditable.
- Each run tests **both** the CRAN release and the GitHub dev head
  (`R CMD INSTALL` of a fresh clone into a separate library; the report
  header prints which library the loaded package actually came from).
- Floating-point comparisons use relative tolerance 1e-12 unless a probe
  documents a wider one at the call site with its justification (e.g.
  fisher.alpha is a numerical root; the solver's own convergence is ~1e-6).
- A probe failure is **not** published as a vegan bug here first: it goes
  to the maintainers privately/upstream before anything else. The suite's
  job is a standing signal, not a gotcha.
- Cadence: weekly, plus a run when a new vegan release appears.

## Probe coverage

**Tier 1 — closed-form definitional oracles (live):**

| Probe | Oracle |
|---|---|
| `diversity()` | Shannon/Simpson/inv-Simpson by hand loop + closed-form edge cases |
| `vegdist()` | Bray, Jaccard (independent Ružička route), Euclidean, Kulczynski, Gower pairwise from textbook formulas; designdist triangulation |
| `decostand()` | total/max/freq/normalize/hellinger/chi.square/log + wisconsin, each as an explicit sweep |
| `rarefy()` | hypergeometric formula **and** brute-force enumeration of all subsamples of a 10-individual community |
| `estimateR()`/`specpool()` | Chao1 (bias-corrected), ACE, Chao2, jackknife 1+2, bootstrap from the documented formulas |
| `renyi()`/`tsallis()` | direct formulas across scales + internal consistency with `diversity()` |
| `fisher.alpha()` | solution-substitution: S = α·log(1 + N/α) |

**Tier 2 — exhaustive-enumeration oracles (planned):** `adonis2`,
`anosim`, `mrpp`, `mantel`, `protest` with *all* permutations enumerated
on tiny designs, exact p-values compared.

**Tier 3 — algebraic identities (planned):** `rda`/`cca` eigenvalues vs
hand SVD + inertia accounting, `betadisper` centroid distances recomputed
from PCoA coordinates, stress of `metaMDS`'s returned configuration,
`specaccum(method="exact")` vs brute force.

Deliberately excluded: null-model internals (`oecosimu`) and anything
where an oracle would have to re-implement vegan's own algorithm —
no independence, no oracle.

## Running it

See the header of [`run-suite.sh`](run-suite.sh) for the exact container
recipe. `Rscript run-probes.R <label>` runs everything against whatever
vegan is first on the library path; exit status 1 on any failure.

## Status

Tier 1 (105 checks) green against CRAN 2.8.0 and dev head as of the
latest report in `reports/`. Maintainers: this suite is strictly opt-in
courtesy infrastructure — if any part of it is unwelcome, say so and it
changes or stops.
