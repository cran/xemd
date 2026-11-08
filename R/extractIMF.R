# First intrinsic mode function E_1 of a signal by sifting (emd.model:9-82,
# :104-226; ceemdan.model:9-15): h_0 is the signal and h_{j+1} = h_j - m_j,
# with m_j the mean of the envelopes of h_j. The extrema and envelopes of each
# state j are computed once and the state is examined in this order: the rule
# accepts it when the rule is defined at j (converged TRUE); j equals fixedSift
# under "fixed" (converged NA); j equals maxSift under an adaptive rule, or the
# envelopes of h_j do not exist (converged FALSE); otherwise one subtraction.
# Returns NULL when h_0 has at most one extremum or no envelopes. sifting is
# validated by the public entry that receives it, not here.
.extractIMF <- function(signal, sifting) {
    Rule <- sifting$stopRule
    Current <- signal
    Extrema <- .findExtrema(Current)
    if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) return(NULL)
    Envelopes <- .buildEnvelopes(Current, extrema = Extrema, boundary = sifting$boundary)
    if (is.null(Envelopes)) return(NULL)
    j <- 0L
    repeat {
        Counts <- c(nrow(Extrema$maxima) + nrow(Extrema$minima), Extrema$crossings)
        Mean <- if (!is.null(Envelopes)) (Envelopes$upper + Envelopes$lower) / 2
        Accepted <- FALSE
        if (Rule == "SD" && j >= 1L) {
            SD <- .measureSD(previous = Previous, current = Current)
            Accepted <- is.finite(SD) && SD < sifting$tolerance
        }
        if (Rule == "Snumber" && j >= 1L) {
            Stable <- if (j >= 2L && identical(Counts, Counts.previous)) Stable + 1L else 1L
            Accepted <- Stable >= sifting$stableSift && abs(Counts[1L] - Counts[2L]) <= 1L
        }
        if (Rule == "doubleThreshold" && !is.null(Envelopes)) {
            Deviation <- .measureDeviation(mean = Mean,
                                           amplitude = (Envelopes$upper - Envelopes$lower) / 2)
            Accepted <- abs(Counts[1L] - Counts[2L]) <= 1L &&
                .meetsDoubleThreshold(Deviation, thresholds = sifting$thresholds,
                                      fraction = sifting$fraction)
        }
        if (Accepted) return(list(mode = Current, iterations = j, converged = TRUE))
        if (Rule == "fixed" && j == sifting$fixedSift) {
            return(list(mode = Current, iterations = j, converged = NA))
        }
        if ((Rule != "fixed" && j == sifting$maxSift) || is.null(Envelopes)) {
            return(list(mode = Current, iterations = j, converged = FALSE,
                        stopReason = if (is.null(Envelopes)) "insufficientExtrema" else "iterationLimit"))
        }
        Previous <- Current
        Counts.previous <- Counts
        Current <- Current - Mean
        j <- j + 1L
        Extrema <- .findExtrema(Current)
        Envelopes <- .buildEnvelopes(Current, extrema = Extrema, boundary = sifting$boundary)
    }
}
