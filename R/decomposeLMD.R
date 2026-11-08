# Local mean decomposition (lmd.model:51-93) by successive extraction of
# product functions: r_0 is the signal; before extraction p the decomposition
# ends with "terminalResidue" when r_{p-1} has at most one extremum, a residue
# constant or without oscillations (lmd.model:93), and with "modeLimit" when
# p - 1 equals maxModes; otherwise PF_p is the product function of r_{p-1}
# (.extractProduct) and r_p = r_{p-1} - PF_p (eq-lmd-residuo). Product
# functions are not intrinsic mode functions. Positions are sample indices,
# and the decomposition runs on the bare samples, so no attribute of signal
# (names, tsp, class) reaches the outputs.
decomposeLMD <- function(signal, maxIterations, envelopeTolerance, changeTolerance,
                         maxSmoothPasses, maxModes = Inf) {
    isCount <- function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
    }
    isTolerance <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
    stopifnot(
        "signal must be a numeric vector of finite values without dim" =
            is.numeric(signal) && is.null(dim(signal)) && all(is.finite(signal)),
        "maxIterations must be a whole number not below one" = isCount(maxIterations),
        "envelopeTolerance must be a finite number not below zero" = isTolerance(envelopeTolerance),
        "changeTolerance must be a finite number not below zero" = isTolerance(changeTolerance),
        "maxSmoothPasses must be a whole number not below one" = isCount(maxSmoothPasses),
        "maxModes must be a whole number not below one or Inf" =
            is.numeric(maxModes) && length(maxModes) == 1L && !is.na(maxModes) &&
            maxModes >= 1 && maxModes == floor(maxModes))
    signal <- as.vector(signal)
    Residue <- signal
    Modes <- list()
    repeat {
        Extrema <- .findExtrema(Residue)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == maxModes) {
            StopReason <- "modeLimit"
            break
        }
        Extraction <- .extractProduct(Residue, maxIterations = maxIterations,
                                      envelopeTolerance = envelopeTolerance,
                                      changeTolerance = changeTolerance,
                                      maxSmoothPasses = maxSmoothPasses)
        Modes[[length(Modes) + 1L]] <- Extraction
        Residue <- Residue - Extraction$mode
    }
    list(
        modes = matrix(as.double(unlist(lapply(Modes, `[[`, "mode"), use.names = FALSE)),
                       nrow = length(signal), ncol = length(Modes)),
        residue = Residue,
        iterations = vapply(Modes, `[[`, integer(1L), "iterations"),
        converged = vapply(Modes, `[[`, logical(1L), "converged"),
        stopReason = StopReason
    )
}
