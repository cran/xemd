# Complete ensemble empirical mode decomposition with adaptive noise
# (ceemdan.model:98-161). The columns w^(m) of the supplied unit-variance
# noise are decomposed once with sifting, eta_k^(m) = E_k(w^(m)), and
# eta_0^(m) = w^(m). From r_0 = signal, before stage k the common residue
# r_{k-1} ends the decomposition with "terminalResidue" when it has at most
# one extremum, with "modeLimit" when k - 1 equals maxModes and with
# "insufficientExtrema" when it has no envelopes. Stage k perturbs r_{k-1}
# with p^(m) = beta_{k-1}^(m) eta_{k-1}^(m), at the constant ratio
# beta_{k-1}^(m) = amplitudeRatio sd(r_{k-1}) / sd(eta_{k-1}^(m)) of each stage
# (ceemdan.m of Colominas for the 2011 formulation; libeemd ceemdan.c), sd of
# divisor N - 1; a realization without eta_{k-1} has a null perturbation and
# stays in the divisor. d_k is the average first mode of r_{k-1} + p^(m)
# (eq-completas-ceemdan-primera, eq-completas-ceemdan-etapa) and
# r_k = r_{k-1} - d_k (eq-completas-ceemdan-recursion). A noise mode of zero
# deviation leaves its scale undefined. amplitudeRatio = 0.2 as in
# decomposeEEMD.
decomposeCEEMDAN <- function(signal, noise, amplitudeRatio = 0.2, sifting) {
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
        Perturbations <- vapply(seq_len(ncol(noise)), function(m) {
            if (k > 1L && ncol(NoiseModes[[m]]) < k - 1L) return(numeric(length(signal)))
            Mode <- if (k == 1L) noise[, m] else NoiseModes[[m]][, k - 1L]
            Deviation <- stats::sd(Mode)
            if (Deviation == 0) {
                stop(errorCondition("a noise mode has zero standard deviation: its scale is undefined",
                                    class = c("xemdZeroEnergy", "moduleFailure")))
            }
            Scale * Mode / Deviation
        }, numeric(length(signal)))
        Average <- .averageFirstMode(Residue, perturbations = Perturbations, sifting = sifting)
        Modes[[k]] <- Average$mode
        Iterations[k] <- Average$iterations
        Converged[k] <- Average$converged
        Residue <- Residue - Average$mode
    }
    list(modes = matrix(as.double(unlist(Modes, use.names = FALSE)), nrow = length(signal), ncol = length(Modes)),
         residue = Residue, iterations = Iterations, converged = Converged, stopReason = StopReason)
}
