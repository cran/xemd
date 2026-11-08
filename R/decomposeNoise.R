# Modes eta_k^(m) = E_k(w^(m)) of each noise realization (ceemdan.model:29,
# :260), decomposed with the configuration of the signal: a list of M
# matrices samples x K_m, in the order of the columns of noise.
.decomposeNoise <- function(noise, sifting) {
    lapply(seq_len(ncol(noise)), function(m) decomposeEMD(noise[, m], sifting = sifting)$modes)
}
