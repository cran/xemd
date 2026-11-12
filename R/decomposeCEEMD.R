# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Complementary ensemble empirical mode decomposition (ceemdan.model:31-59).
# For each column z^(m) of the supplied unit-variance noise, the inputs
# x + a z^(m) and x - a z^(m), with a = amplitudeRatio * sd(x) and sd of
# divisor N - 1 (eq-completas-par), are decomposed by scalar sifting, the
# positive input first. Each stops before an incomplete extraction and keeps
# its input in the residue. The 2M decompositions are averaged
# by extraction order with zero extension and the divisor 2M of
# eq-completas-modo-complementario, and residue is their mean residue, so the
# perturbations cancel in the reconstruction (eq-completas-cancelacion-par).
# amplitudeRatio = 0.2 is the ratio common to libeemd and emd 0.8.1, as in
# decomposeEEMD.
decomposeCEEMD <- function(signal, noise, amplitudeRatio = 0.2, sifting) {
    stopifnot(
        "signal must be a numeric vector of at least two finite values without dim" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 2L && all(is.finite(signal)),
        "noise must be a finite numeric matrix with one row per sample of signal" =
            is.matrix(noise) && is.numeric(noise) && nrow(noise) == length(signal) &&
            ncol(noise) >= 1L && all(is.finite(noise)),
        "amplitudeRatio must be a finite number not below zero" =
            is.numeric(amplitudeRatio) && length(amplitudeRatio) == 1L &&
            is.finite(amplitudeRatio) && amplitudeRatio >= 0)
    signal <- as.vector(signal)
    noise <- matrix(as.double(noise), nrow = nrow(noise))
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule == "consecutive") {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    Perturbations <- amplitudeRatio * stats::sd(signal) * noise
    Decompositions <- unlist(lapply(seq_len(ncol(Perturbations)), function(m) {
        list(.decomposeScalar(signal + Perturbations[, m], sifting = sifting, acceptIncomplete = FALSE),
             .decomposeScalar(signal - Perturbations[, m], sifting = sifting, acceptIncomplete = FALSE))
    }), recursive = FALSE)
    OUT <- .averageDecompositions(.alignRealizations(Decompositions), divisor = "realizations")
    OUT$termination <- do.call(rbind, lapply(seq_along(Decompositions), function(m) {
        data.frame(realization = (m + 1L) %/% 2L, sign = if (m %% 2L) 1L else -1L,
                   Decompositions[[m]]$termination)
    }))
    OUT
}
