# decostand(): every method checked here has a one-line definition in the
# documentation; each is recomputed with explicit sweeps in base R.
# wisconsin() = species standardized by their maxima, then sites by totals.

comm <- as.matrix(read.csv("data/community.csv"))

total_hand     <- comm / rowSums(comm)
max_hand       <- sweep(comm, 2, apply(comm, 2, max), "/")
freq_hand      <- sweep(sweep(comm, 2, colSums(comm), "/"), 2, colSums(comm > 0), "*")
normalize_hand <- comm / sqrt(rowSums(comm^2))
hellinger_hand <- sqrt(comm / rowSums(comm))
chi_hand       <- comm * sqrt(sum(comm)) / outer(rowSums(comm), sqrt(colSums(comm)))
log_hand       <- ifelse(comm > 0, log(comm, base = 2) + 1, 0)
wis_step1      <- sweep(comm, 2, apply(comm, 2, max), "/")
wisconsin_hand <- wis_step1 / rowSums(wis_step1)

check_close("decostand total",     total_hand,     decostand(comm, "total"))
check_close("decostand max",       max_hand,       decostand(comm, "max"))
check_close("decostand freq",      freq_hand,      decostand(comm, "freq"))
check_close("decostand normalize", normalize_hand, decostand(comm, "normalize"))
check_close("decostand hellinger", hellinger_hand, decostand(comm, "hellinger"))
check_close("decostand chi.square", chi_hand,      decostand(comm, "chi.square"))
check_close("decostand log base 2", log_hand,      decostand(comm, "log"))
check_close("wisconsin double standardization", wisconsin_hand, wisconsin(comm))
