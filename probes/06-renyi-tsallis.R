# renyi()/tsallis(): direct formulas across scales, plus internal
# consistency with diversity() at the coinciding orders (scale 1 = Shannon
# for both; renyi scale 0 = log richness; hill numbers = exp(renyi)).

x <- c(12, 5, 1, 2, 8, 4)
p <- x / sum(x)
scales <- c(0, 0.5, 2, 4)

ren_hand <- sapply(scales, function(a) log(sum(p^a)) / (1 - a))
check_close("renyi formula, scales 0/0.5/2/4", ren_hand,
            as.vector(as.matrix(renyi(x, scales = scales))))
check_close("renyi scale 1 = shannon", unname(diversity(x, "shannon")),
            as.vector(as.matrix(renyi(x, scales = 1))))
check_close("renyi hill = exp(renyi)", exp(ren_hand),
            as.vector(as.matrix(renyi(x, scales = scales, hill = TRUE))))

tsa_hand <- sapply(scales, function(a) (1 - sum(p^a)) / (a - 1))
check_close("tsallis formula, scales 0/0.5/2/4", tsa_hand,
            as.vector(as.matrix(tsallis(x, scales = scales))))
check_close("tsallis scale 1 = shannon", unname(diversity(x, "shannon")),
            as.vector(as.matrix(tsallis(x, scales = 1))))
check_close("tsallis scale 0 = S - 1", sum(x > 0) - 1,
            as.vector(as.matrix(tsallis(x, scales = 0))))
