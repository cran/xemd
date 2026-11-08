# Consecutive post-processing of ensemble modes D_1, ..., D_K, K >= 2
# (eemd.model:231-265). With u_1 = D_1 + D_2, stage i extracts
# a_i = E_1(u_i) (eq-eemd-posprocesado-extraccion), keeps the remainder
# q_i = u_i - a_i (eq-eemd-posprocesado-remanente) and adds the next mode,
# u_{i+1} = q_i + D_{i+2} (eq-eemd-posprocesado-cadena). The completed chain
# returns K - 1 modes and the residue q_{K-1} with "chainComplete", and
# sum_i a_i + q_{K-1} = sum_j D_j (eq-eemd-posprocesado-reconstruccion). When
# E_1(u_i) does not exist the chain ends with a_1, ..., a_{i-1} and the
# residue u_i + sum_{j > i + 1} D_j, which keeps the sum, with
# "insufficientExtrema". The chain has K - 1 extractions, so maxModes of
# sifting does not apply.
resiftModes <- function(modes, sifting) {
    stopifnot("modes must be a finite numeric matrix with at least two columns" =
                  is.matrix(modes) && is.numeric(modes) && ncol(modes) >= 2L && all(is.finite(modes)))
    modes <- matrix(as.double(modes), nrow = nrow(modes))
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule == "consecutive") {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    K <- ncol(modes)
    Extractions <- list()
    StopReason <- "chainComplete"
    Chain <- modes[, 1L] + modes[, 2L]
    for (i in seq_len(K - 1L)) {
        Extraction <- .extractIMF(Chain, sifting = sifting)
        if (is.null(Extraction)) {
            StopReason <- "insufficientExtrema"
            Chain <- Chain + rowSums(modes[, -seq_len(i + 1L), drop = FALSE])
            break
        }
        Extractions[[i]] <- Extraction
        Chain <- Chain - Extraction$mode
        if (i < K - 1L) Chain <- Chain + modes[, i + 2L]
    }
    list(modes = matrix(as.double(unlist(lapply(Extractions, `[[`, "mode"), use.names = FALSE)),
                        nrow = nrow(modes), ncol = length(Extractions)),
         residue = Chain,
         iterations = vapply(Extractions, `[[`, integer(1L), "iterations"),
         converged = vapply(Extractions, `[[`, logical(1L), "converged"),
         stopReason = StopReason)
}
