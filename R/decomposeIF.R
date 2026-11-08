# Iterative filtering (lmd.model:224-299) by successive extraction: r_0 is
# the signal; before extraction k the decomposition ends with
# "terminalResidue" when r_{k-1} has at most one extremum and with
# "modeLimit" when k - 1 equals maxModes. Otherwise the half width l_k follows
# from the extrema of r_{k-1} (.selectHalfWidth) and stays fixed during the
# extraction, the supplied prototype filter is resampled to 2 l_k + 1 weights
# (.resampleFilter), c_k = f^(Q_k) of the iteration with that filter
# (eq-if-iteracion) and r_k = r_{k-1} - c_k. boundary "zero" extends the
# record with zeros; "periodic" extends it periodically and evaluates the
# iteration spectrally, which is FIF (lmd.model:428-442). tolerance bounds the
# relative change of eq-if-parada-relativa between consecutive iterates, a
# ratio of norms. Positions are sample indices, and the decomposition runs on
# the bare samples, so no attribute of signal (names, tsp, class) reaches the
# outputs. The components are not shown to satisfy the counting conditions of
# intrinsic mode functions (lmd.model:299).
decomposeIF <- function(signal, filter, boundary, tolerance, maxIterations, stretch,
                        maxModes = Inf) {
    stopifnot(
        "signal must be a numeric vector of finite values without dim" =
            is.numeric(signal) && is.null(dim(signal)) && all(is.finite(signal)),
        "filter must be a numeric vector of finite values of odd length, at least 3" =
            is.numeric(filter) && is.null(dim(filter)) && length(filter) >= 3L &&
            length(filter) %% 2L == 1L && all(is.finite(filter)),
        "boundary must be \"periodic\" or \"zero\"" =
            is.character(boundary) && length(boundary) == 1L && boundary %in% c("periodic", "zero"),
        "tolerance must be a finite number greater than zero and less than one" =
            is.numeric(tolerance) && length(tolerance) == 1L && is.finite(tolerance) &&
            tolerance > 0 && tolerance < 1,
        "maxIterations must be a whole number not below one" =
            is.numeric(maxIterations) && length(maxIterations) == 1L && is.finite(maxIterations) &&
            maxIterations >= 1 && maxIterations == floor(maxIterations),
        "stretch must be a finite number greater than zero" =
            is.numeric(stretch) && length(stretch) == 1L && is.finite(stretch) && stretch > 0,
        "maxModes must be a whole number not below one or Inf" =
            is.numeric(maxModes) && length(maxModes) == 1L && !is.na(maxModes) &&
            maxModes >= 1 && maxModes == floor(maxModes))
    signal <- as.vector(signal)
    iterate <- switch(boundary, periodic = .iterateFluctuationPeriodic, zero = .iterateFluctuationZero)
    Residue <- signal
    Modes <- list()
    HalfWidth <- NULL
    repeat {
        Extrema <- .findExtrema(Residue)
        ExtremaCount <- nrow(Extrema$maxima) + nrow(Extrema$minima)
        if (ExtremaCount <= 1L) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == maxModes) {
            StopReason <- "modeLimit"
            break
        }
        HalfWidth <- .selectHalfWidth(extremaCount = ExtremaCount, sampleCount = length(signal),
                                      stretch = stretch, previous = HalfWidth)
        Extraction <- iterate(Residue, filter = .resampleFilter(filter, halfWidth = HalfWidth),
                              tolerance = tolerance, maxIterations = maxIterations)
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
