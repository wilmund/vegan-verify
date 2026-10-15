# estimateR() / specpool(): richness extrapolators recomputed from the
# formulas vegan's own documentation states (?specpool Details) — the
# documentation is the independent definition here. Data constructed so
# that singleton/doubleton counts are visible by eye.

## estimateR (abundance-based, single site):
## bias-corrected Chao1 = S0 + a1*(a1-1)/(2*(a2+1)), no small-sample corr;
## ACE with rare-species threshold 10.
x <- c(1, 1, 1, 2, 2, 3, 4, 12, 15)   # S0 = 9, a1 = 3, a2 = 2
est <- estimateR(x)

S0 <- sum(x > 0)
a  <- tabulate(x)                      # a[i] = number of species with abundance i
chao1_hand <- S0 + a[1] * (a[1] - 1) / (2 * (a[2] + 1))

rare    <- x[x <= 10]
S_rare  <- length(rare); S_abund <- sum(x > 10); N_rare <- sum(rare)
C_ace   <- 1 - a[1] / N_rare
i       <- seq_len(min(10, length(a)))
gamma2  <- max(S_rare / C_ace * sum(i * (i - 1) * a[i]) / (N_rare * (N_rare - 1)) - 1, 0)
ace_hand <- S_abund + S_rare / C_ace + a[1] / C_ace * gamma2

check_close("estimateR S.obs",   S0,         unname(est["S.obs"]))
check_close("estimateR S.chao1", chao1_hand, unname(est["S.chao1"]))
check_close("estimateR S.ACE",   ace_hand,   unname(est["S.ACE"]))

# a2 = 0 case: bias-corrected form must still apply (doc: no basic-form switch
# for estimateR).
x2 <- c(1, 1, 5, 8)                    # S0 = 4, a1 = 2, a2 = 0
check_close("estimateR S.chao1 with a2=0",
            4 + 2 * 1 / 2, unname(estimateR(x2)["S.chao1"]))

## specpool (incidence-based, smallsample = TRUE default):
## chao  = S0 + a1^2/(2*a2) * (N-1)/N        (a2 > 0 here)
## jack1 = S0 + a1*(N-1)/N
## jack2 = S0 + a1*(2N-3)/N - a2*(N-2)^2/N/(N-1)
## boot  = S0 + sum((1-p_i)^N), p_i = occurrence frequency of species i
inc <- matrix(c(1, 1, 1, 1,
                1, 0, 0, 0,
                0, 1, 0, 0,
                1, 1, 0, 0,
                0, 0, 1, 1,
                1, 1, 1, 0,
                0, 1, 1, 1), nrow = 4,
              dimnames = list(paste0("site", 1:4), paste0("sp", 1:7)))
N  <- nrow(inc)
occ <- colSums(inc > 0)
S0 <- sum(occ > 0); a1 <- sum(occ == 1); a2 <- sum(occ == 2)
p  <- occ / N
pool <- specpool(inc)

check_close("specpool chao",  S0 + a1^2 / (2 * a2) * (N - 1) / N, pool$chao)
check_close("specpool jack1", S0 + a1 * (N - 1) / N,              pool$jack1)
check_close("specpool jack2", S0 + a1 * (2 * N - 3) / N - a2 * (N - 2)^2 / N / (N - 1),
            pool$jack2)
check_close("specpool boot",  S0 + sum((1 - p)^N),                pool$boot)
