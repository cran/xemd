# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Ensemble empirical mode decomposition (eemd.model:17-79, :185-207). Each
# realization x + w^(m), with w^(m) = amplitudeRatio * sd(x) * z^(m), z^(m)
# the column m of the supplied unit-variance noise and sd of divisor N - 1
# (eq-eemd-entrada, eq-eemd-amplitud-relativa), is decomposed by
# scalar sifting (eq-eemd-realizacion), stopping before an incomplete
# extraction and retaining its input in the residue, in the order of the
# columns. The modes are averaged by extraction order with zero extension and
# the divisor M for "realizations" (eq-eemd-divisores) or the divisor M_j of
# the realizations that have the mode for "available"
# (eq-eemd-media-disponibles); residue is the mean of the residues
# (eq-eemd-media-residuo) and meanNoise the mean of the added perturbations
# (eq-eemd-ruido-medio). A realization without modes adds zeros to the modes
# and its input to the residue. amplitudeRatio = 0.2 is the ratio to the
# standard deviation of the signal common to libeemd and emd 0.8.1
# (Rlibeemd/R/eemd.R:56; emd/sift.py:693, :803).
decomposeEEMD <- function(signal, noise, amplitudeRatio = 0.2, divisor, sifting) {
    stopifnot(
        "signal must be a numeric vector of at least two finite values without dim" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 2L && all(is.finite(signal)),
        "noise must be a finite numeric matrix with one row per sample of signal" =
            is.matrix(noise) && is.numeric(noise) && nrow(noise) == length(signal) &&
            ncol(noise) >= 1L && all(is.finite(noise)),
        "amplitudeRatio must be a finite number not below zero" =
            is.numeric(amplitudeRatio) && length(amplitudeRatio) == 1L &&
            is.finite(amplitudeRatio) && amplitudeRatio >= 0,
        "divisor must be realizations or available" =
            is.character(divisor) && length(divisor) == 1L && divisor %in% c("realizations", "available"))
    signal <- as.vector(signal)
    noise <- matrix(as.double(noise), nrow = nrow(noise))
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule == "consecutive") {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    Perturbations <- amplitudeRatio * stats::sd(signal) * noise
    Decompositions <- lapply(seq_len(ncol(Perturbations)), function(m) {
        .decomposeScalar(signal + Perturbations[, m], sifting = sifting, acceptIncomplete = FALSE)
    })
    OUT <- .averageDecompositions(.alignRealizations(Decompositions), divisor = divisor)
    list(modes = OUT$modes, residue = OUT$residue, meanNoise = rowMeans(Perturbations),
         iterations = OUT$iterations, converged = OUT$converged,
         termination = do.call(rbind, lapply(seq_along(Decompositions), function(m) {
             data.frame(realization = m, Decompositions[[m]]$termination)
         })))
}
