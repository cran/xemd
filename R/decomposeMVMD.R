# Multivariate variational mode decomposition (vmd.extensiones.es.qmd:1-153)
# of the N x C matrix signal, sampled every samplingInterval seconds, into
# length(initialCenters) modes whose C coordinates share one center. Frequencies
# are in Hz at the interface and in rad/s inside, where alpha has units of s^2
# (eq-vmd-escala-alpha). tau > 0 solves the constrained problem and tau = 0 the
# penalized one with zero multipliers (vmd.extensiones.es.qmd:143);
# fixedZeroCenter keeps the first center at 0 while its coordinates are still
# updated, as DC = true in the authors' MVMD.m 1.0.1. Sweeps start from zero
# modes and multipliers and stop when the relative change over modes and
# channels (eq-mvmd-parada) falls below tolerance, which cannot happen before
# the second sweep, or after maxIterations sweeps; the default 499 is the limit
# of MVMD.m 1.0.1 (500 states).
decomposeMVMD <- function(signal, samplingInterval, initialCenters, alpha, tau,
                          fixedZeroCenter, tolerance, maxIterations = 499) {
    stopifnot(
        "signal must be a numeric matrix of finite values with at least one row and one column" =
            is.numeric(signal) && is.matrix(signal) && nrow(signal) >= 1L &&
            ncol(signal) >= 1L && all(is.finite(signal)),
        "samplingInterval must be a finite number greater than zero" =
            is.numeric(samplingInterval) && length(samplingInterval) == 1L &&
            is.finite(samplingInterval) && samplingInterval > 0,
        "initialCenters must be finite frequencies from 0 to 1 / (2 * samplingInterval)" =
            is.numeric(initialCenters) && is.null(dim(initialCenters)) &&
            length(initialCenters) >= 1L && all(is.finite(initialCenters)) &&
            all(initialCenters >= 0) &&
            all(initialCenters <= 1 / (2 * samplingInterval)),
        "alpha must be a finite number greater than zero" =
            is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
            alpha > 0,
        "tau must be a finite number not below zero" =
            is.numeric(tau) && length(tau) == 1L && is.finite(tau) && tau >= 0,
        "fixedZeroCenter must be TRUE or FALSE" =
            is.logical(fixedZeroCenter) && length(fixedZeroCenter) == 1L &&
            !is.na(fixedZeroCenter),
        "fixedZeroCenter requires initialCenters[1] == 0" =
            !fixedZeroCenter || initialCenters[1L] == 0,
        "tolerance must be a finite number greater than zero" =
            is.numeric(tolerance) && length(tolerance) == 1L &&
            is.finite(tolerance) && tolerance > 0,
        "maxIterations must be a whole number not below one" =
            is.numeric(maxIterations) && length(maxIterations) == 1L &&
            is.finite(maxIterations) && maxIterations >= 1 &&
            maxIterations == floor(maxIterations))
    # The outputs inherit no attribute of signal other than its dimensions.
    signal <- matrix(as.vector(signal), nrow = nrow(signal))
    Transform <- .transformReflection(signal = signal,
                                      samplingInterval = samplingInterval)
    BinCount <- nrow(Transform$spectrum)
    ModeCount <- length(initialCenters)
    State <- list(modeSpectra = array(0i, dim = c(BinCount, ModeCount,
                                                  ncol(signal))),
                  angularCenters = 2 * pi * initialCenters,
                  multiplier = matrix(0i, nrow = BinCount, ncol = ncol(signal)))
    for (Iteration in seq_len(maxIterations)) {
        Previous <- State$modeSpectra
        State <- .sweepVMD(spectrum = Transform$spectrum,
                           modeSpectra = Previous,
                           angularCenters = State$angularCenters,
                           multiplier = State$multiplier,
                           angularFrequency = Transform$angularFrequency,
                           alpha = alpha, tau = tau,
                           fixedZeroCenter = fixedZeroCenter)
        Converged <- .measureSweepChange(current = State$modeSpectra,
                                         previous = Previous) < tolerance
        if (Converged) break
    }
    # The T x (K C) matrix of the mode spectra restores to samples x modes x
    # channels.
    Modes <- array(.restoreObserved(
        spectrum = matrix(State$modeSpectra, nrow = BinCount),
        observed = Transform$observed), dim = c(nrow(signal), ModeCount,
                                                ncol(signal)))
    list(modes = Modes,
         fitResidual = signal - rowSums(aperm(Modes, c(1L, 3L, 2L)),
                                        dims = 2L),
         centers = State$angularCenters / (2 * pi),
         iterations = rep(Iteration, times = ModeCount),
         converged = rep(Converged, times = ModeCount),
         stopReason = if (Converged) "tolerance" else "iterationLimit")
}
