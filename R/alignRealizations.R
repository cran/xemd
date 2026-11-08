# Zero extension of M decompositions to K_max modes (eq-eemd-divisores,
# eemd.model:195-199): the modes of realization m fill slice m of an array
# samples x K_max x M, with zero modes where it has fewer than K_max; the
# residues are the columns of a matrix samples x M; counts and flags are
# matrices K_max x M with zero counts and NA flags for absent modes; modeCount
# holds K_m, from which M_j = #{m : K_m >= j} follows.
.alignRealizations <- function(decompositions) {
    Counts <- vapply(decompositions, function(x) ncol(x$modes), integer(1L))
    n <- length(decompositions[[1L]]$residue)
    M <- length(decompositions)
    K <- max(Counts)
    OUT <- list(modes = array(0, dim = c(n, K, M)),
                residue = matrix(vapply(decompositions, `[[`, numeric(n), "residue"), nrow = n, ncol = M),
                iterations = matrix(0L, nrow = K, ncol = M),
                converged = matrix(NA, nrow = K, ncol = M),
                modeCount = Counts)
    for (m in seq_len(M)) {
        IDX <- seq_len(Counts[m])
        OUT$modes[, IDX, m] <- decompositions[[m]]$modes
        OUT$iterations[IDX, m] <- decompositions[[m]]$iterations
        OUT$converged[IDX, m] <- decompositions[[m]]$converged
    }
    OUT
}
