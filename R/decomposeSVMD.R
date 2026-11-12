# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Successive variational mode decomposition (vmd.extensiones.es.qmd:165-350) of
# a signal sampled every samplingInterval seconds. Mode L is extracted from the
# whole signal with the centers of modes 1, ..., L - 1 fixed
# (eq-svmd-particion), starting from initialCenter in Hz; frequencies are in Hz
# at the interface and in rad/s inside, where alpha and maxAlpha have units of
# s^2. Modes keep the extraction order. stopping ends the extraction: rule
# "modeCount" after stopping$modeCount modes; rule "noisePower" when the mean
# square of the residual on the observed samples (eq-svmd-potencia-residual)
# meets eq-svmd-parada-exterior with stopping$noisePower in squared signal
# units and stopping$noiseTolerance, or after stopping$maxModes modes (default
# Inf). The defaults initialCenter = 0 and maxIterations = 299 are those of the
# authors' svmd.m 1.1.1.
decomposeSVMD <- function(signal, samplingInterval, initialCenter = 0, alpha,
                          maxAlpha, tau, tolerance, maxIterations = 299,
                          stopping) {
    stopifnot(
        "signal must be a non-empty numeric vector of finite values" =
            is.numeric(signal) && is.null(dim(signal)) &&
            length(signal) >= 1L && all(is.finite(signal)),
        "samplingInterval must be a finite number greater than zero" =
            is.numeric(samplingInterval) && length(samplingInterval) == 1L &&
            is.finite(samplingInterval) && samplingInterval > 0,
        "initialCenter must be a finite frequency from 0 to 1 / (2 * samplingInterval)" =
            is.numeric(initialCenter) && length(initialCenter) == 1L &&
            is.finite(initialCenter) && initialCenter >= 0 &&
            initialCenter <= 1 / (2 * samplingInterval),
        "alpha must be a finite number greater than zero" =
            is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
            alpha > 0,
        "maxAlpha must be a finite number not below alpha" =
            is.numeric(maxAlpha) && length(maxAlpha) == 1L &&
            is.finite(maxAlpha) && maxAlpha >= alpha,
        "tau must be a finite number not below zero" =
            is.numeric(tau) && length(tau) == 1L && is.finite(tau) && tau >= 0,
        "tolerance must be a finite number greater than zero" =
            is.numeric(tolerance) && length(tolerance) == 1L &&
            is.finite(tolerance) && tolerance > 0,
        "maxIterations must be a whole number not below one" =
            is.numeric(maxIterations) && length(maxIterations) == 1L &&
            is.finite(maxIterations) && maxIterations >= 1 &&
            maxIterations == floor(maxIterations),
        "stopping$stopRule must be \"noisePower\" or \"modeCount\"" =
            is.list(stopping) && is.character(stopping[["stopRule"]]) &&
            length(stopping[["stopRule"]]) == 1L &&
            stopping[["stopRule"]] %in% c("noisePower", "modeCount"))
    if (stopping[["stopRule"]] == "noisePower" &&
        !"maxModes" %in% names(stopping)) {
        stopping$maxModes <- Inf
    }
    Fields <- list(noisePower = c("stopRule", "noisePower", "noiseTolerance",
                                  "maxModes"),
                   modeCount = c("stopRule", "modeCount"))
    Fields <- Fields[[stopping[["stopRule"]]]]
    stopifnot("stopping must hold exactly the fields of its stopRule" =
                  length(stopping) == length(Fields) &&
                  setequal(names(stopping), Fields))
    if (stopping[["stopRule"]] == "modeCount") {
        stopifnot("stopping$modeCount must be a whole number not below one" =
                      is.numeric(stopping[["modeCount"]]) &&
                      length(stopping[["modeCount"]]) == 1L &&
                      is.finite(stopping[["modeCount"]]) &&
                      stopping[["modeCount"]] >= 1 &&
                      stopping[["modeCount"]] == floor(stopping[["modeCount"]]))
    }
    if (stopping[["stopRule"]] == "noisePower") {
        stopifnot(
            "stopping$noisePower must be a finite number greater than zero" =
                is.numeric(stopping[["noisePower"]]) &&
                length(stopping[["noisePower"]]) == 1L &&
                is.finite(stopping[["noisePower"]]) &&
                stopping[["noisePower"]] > 0,
            "stopping$noiseTolerance must be a finite number greater than zero" =
                is.numeric(stopping[["noiseTolerance"]]) &&
                length(stopping[["noiseTolerance"]]) == 1L &&
                is.finite(stopping[["noiseTolerance"]]) &&
                stopping[["noiseTolerance"]] > 0,
            "stopping$maxModes must be a whole number not below one or Inf" =
                is.numeric(stopping[["maxModes"]]) &&
                length(stopping[["maxModes"]]) == 1L &&
                !is.na(stopping[["maxModes"]]) &&
                stopping[["maxModes"]] >= 1 &&
                stopping[["maxModes"]] == floor(stopping[["maxModes"]]))
    }
    # The outputs inherit no attribute of signal (names, tsp, class).
    signal <- as.vector(signal)
    Transform <- .transformReflection(signal = signal,
                                      samplingInterval = samplingInterval)
    Modes <- matrix(0, nrow = length(signal), ncol = 0)
    Centers <- numeric()
    Iterations <- integer()
    Converged <- logical()
    StopReason <- NULL
    while (is.null(StopReason)) {
        Mode <- .extractVariationalMode(
            spectrum = Transform$spectrum[, 1L],
            angularFrequency = Transform$angularFrequency,
            previousCenters = Centers, initialCenter = 2 * pi * initialCenter,
            alpha = alpha, maxAlpha = maxAlpha, tau = tau,
            tolerance = tolerance, maxIterations = maxIterations)
        Modes <- cbind(Modes, .restoreObserved(
            spectrum = matrix(Mode$modeSpectrum),
            observed = Transform$observed))
        Centers <- c(Centers, Mode$center)
        Iterations <- c(Iterations, Mode$iterations)
        Converged <- c(Converged, Mode$converged)
        if (stopping[["stopRule"]] == "modeCount" &&
            ncol(Modes) == stopping[["modeCount"]]) {
            StopReason <- "modeCount"
        }
        if (stopping[["stopRule"]] == "noisePower") {
            Power <- mean((signal - rowSums(Modes))^2)
            Met <- abs(stopping[["noisePower"]] - Power) /
                stopping[["noisePower"]] < stopping[["noiseTolerance"]]
            if (Met) StopReason <- "noisePower"
            if (!Met && ncol(Modes) == stopping[["maxModes"]]) {
                StopReason <- "modeLimit"
            }
        }
    }
    list(modes = Modes, fitResidual = signal - rowSums(Modes),
         centers = Centers / (2 * pi), iterations = Iterations,
         converged = Converged, stopReason = StopReason)
}
