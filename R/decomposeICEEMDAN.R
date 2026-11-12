# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Improved complete ensemble empirical mode decomposition with adaptive noise
# (ceemdan.model:7-27, :165-256). The columns w^(m) of the supplied
# unit-variance noise are decomposed once with sifting, eta_k^(m) = E_k(w^(m)).
# From r_0 = signal, before stage k the common residue r_{k-1} ends the
# decomposition with "terminalResidue" when it has at most one extremum,
# with "modeLimit" when k - 1 equals maxModes and with "insufficientExtrema"
# when it has no envelopes. Stage k perturbs r_{k-1} with
# p^(m) = beta_{k-1}^(m) eta_k^(m): beta_0^(m) = amplitudeRatio sd(x) /
# sd(eta_1^(m)) (eq-completas-normalizacion-inicial) and, for k >= 2,
# beta_{k-1} = amplitudeRatio sd(r_{k-1}), divided by sd(eta_k^(m)) when
# normalizeNoise (eq-completas-normalizacion-posterior), sd of divisor N - 1.
# A realization without eta_k has a null perturbation and stays in the
# divisor. With A the average first mode of r_{k-1} + p^(m) and p the mean
# perturbation, r_k = r_{k-1} + p - A is the mean local mean
# (eq-completas-iceemdan-residuo, eq-completas-diferencia-promedios) and
# d_k = r_{k-1} - r_k (eq-completas-iceemdan-modo). A noise mode of zero
# deviation leaves its scale undefined (ceemdan.model:237). The defaults
# follow decomposeEEMD; normalizeNoise has none.
decomposeICEEMDAN <- function(signal, noise, amplitudeRatio = 0.2, normalizeNoise, sifting) {
    stopifnot(
        "signal must be a numeric vector of at least two finite values without dim" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 2L && all(is.finite(signal)),
        "noise must be a finite numeric matrix with one row per sample of signal" =
            is.matrix(noise) && is.numeric(noise) && nrow(noise) == length(signal) &&
            ncol(noise) >= 1L && all(is.finite(noise)),
        "amplitudeRatio must be a finite number not below zero" =
            is.numeric(amplitudeRatio) && length(amplitudeRatio) == 1L &&
            is.finite(amplitudeRatio) && amplitudeRatio >= 0,
        "normalizeNoise must be TRUE or FALSE" =
            is.logical(normalizeNoise) && length(normalizeNoise) == 1L && !is.na(normalizeNoise))
    signal <- as.vector(signal)
    noise <- matrix(as.double(noise), nrow = nrow(noise))
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule == "consecutive") {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    NoiseModes <- .decomposeNoise(noise, sifting = sifting)
    Residue <- signal
    Modes <- list()
    Iterations <- integer()
    Converged <- logical()
    repeat {
        Extrema <- .findExtrema(Residue)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == sifting$maxModes) {
            StopReason <- "modeLimit"
            break
        }
        if (is.null(.buildEnvelopes(Residue, extrema = Extrema, boundary = sifting$boundary))) {
            StopReason <- "insufficientExtrema"
            break
        }
        k <- length(Modes) + 1L
        Scale <- amplitudeRatio * stats::sd(Residue)
        Perturbations <- vapply(NoiseModes, function(x) {
            if (ncol(x) < k) return(numeric(length(signal)))
            if (k > 1L && !normalizeNoise) return(Scale * x[, k])
            Deviation <- stats::sd(x[, k])
            if (Deviation == 0) {
                stop(errorCondition("a noise mode has zero standard deviation: its scale is undefined",
                                    class = c("xemdZeroEnergy", "moduleFailure")))
            }
            Scale * x[, k] / Deviation
        }, numeric(length(signal)))
        Average <- .averageFirstMode(Residue, perturbations = Perturbations, sifting = sifting)
        Next <- Residue + rowMeans(Perturbations) - Average$mode
        Modes[[k]] <- Residue - Next
        Iterations[k] <- Average$iterations
        Converged[k] <- Average$converged
        Residue <- Next
    }
    list(modes = matrix(as.double(unlist(Modes, use.names = FALSE)), nrow = length(signal), ncol = length(Modes)),
         residue = Residue, iterations = Iterations, converged = Converged, stopReason = StopReason)
}
