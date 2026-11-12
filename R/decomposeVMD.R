# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Variational mode decomposition (vmd.model.es.qmd) of a signal sampled every
# samplingInterval seconds into length(initialCenters) modes. Frequencies are
# in Hz at the interface and in rad/s inside, where alpha has units of s^2
# (eq-vmd-escala-alpha). tau > 0 solves the constrained problem and tau = 0
# the penalized one with a zero multiplier (vmd.model.es.qmd:233);
# fixedZeroCenter keeps the first center at 0 (vmd.model.es.qmd:271). Sweeps
# start from zero modes and multiplier and stop when the relative change
# (eq-vmd-parada) falls below tolerance, which cannot happen before the
# second sweep, or after maxIterations sweeps.
decomposeVMD <- function(signal, samplingInterval, initialCenters, alpha, tau,
                         fixedZeroCenter, tolerance, maxIterations) {
    stopifnot(
        "signal must be a non-empty numeric vector of finite values" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 1L &&
            all(is.finite(signal)),
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
    # The outputs inherit no attribute of signal (names, tsp, class).
    signal <- as.vector(signal)
    Transform <- .transformReflection(signal = signal,
                                      samplingInterval = samplingInterval)
    BinCount <- nrow(Transform$spectrum)
    ModeCount <- length(initialCenters)
    State <- list(modeSpectra = array(0i, dim = c(BinCount, ModeCount, 1L)),
                  angularCenters = 2 * pi * initialCenters,
                  multiplier = matrix(0i, nrow = BinCount, ncol = 1L))
    for (Iteration in seq_len(maxIterations)) {
        Previous <- State$modeSpectra
        State <- .sweepVMD(spectrum = Transform$spectrum, modeSpectra = Previous,
                           angularCenters = State$angularCenters,
                           multiplier = State$multiplier,
                           angularFrequency = Transform$angularFrequency,
                           alpha = alpha, tau = tau,
                           fixedZeroCenter = fixedZeroCenter)
        Converged <- .measureSweepChange(current = State$modeSpectra,
                                         previous = Previous) < tolerance
        if (Converged) break
    }
    # With one channel the mode spectra form the T x K matrix of the modes.
    Modes <- .restoreObserved(
        spectrum = matrix(State$modeSpectra, nrow = BinCount),
        observed = Transform$observed)
    list(modes = Modes, fitResidual = signal - rowSums(Modes),
         centers = State$angularCenters / (2 * pi),
         iterations = rep(Iteration, times = ModeCount),
         converged = rep(Converged, times = ModeCount),
         stopReason = if (Converged) "tolerance" else "iterationLimit")
}
