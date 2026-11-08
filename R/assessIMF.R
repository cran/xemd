# Intrinsic mode function conditions of each column (imf.model:5-41): count
# is |N_e - N_z| <= 1 (eq-imf-conteo) with the extrema and crossings of
# .findExtrema; deviation is the double threshold on sigma = |m_e / a_e|
# (eq-imf-desviacion; emd.model:177-208), with m_e and a_e the half sum and
# half difference of the envelopes (eq-imf-media,
# eq-imf-amplitud-envolventes), the same test that accepts a mode of the
# double-threshold sifting. deviation is NA where the boundary rule leaves the
# column without envelopes. The defaults are those of emd.m, the code of the
# authors of the double threshold, shared with sifting
# (emd.rilling.implementation:5, :57).
assessIMF <- function(modes, thresholds = c(0.05, 0.5), fraction = 0.05, boundary = "reflection") {
    if (!is.matrix(modes) || !is.numeric(modes) || !all(is.finite(modes))) {
        stop("'modes' must be a finite numeric matrix")
    }
    if (!is.numeric(thresholds) || length(thresholds) != 2L || !all(is.finite(thresholds)) ||
        thresholds[1L] <= 0 || thresholds[1L] >= thresholds[2L]) {
        stop("'thresholds' must be two finite numbers with 0 < thresholds[1] < thresholds[2]")
    }
    if (!is.numeric(fraction) || length(fraction) != 1L || !is.finite(fraction) ||
        fraction < 0 || fraction >= 1) {
        stop("'fraction' must be a finite number in [0, 1)")
    }
    if (!identical(boundary, "reflection") && !identical(boundary, "extrapolation")) {
        stop("'boundary' must be \"reflection\" or \"extrapolation\"")
    }
    attributes(modes) <- list(dim = dim(modes))
    OUT <- list(count = logical(ncol(modes)), deviation = rep(NA, ncol(modes)))
    for (k in seq_len(ncol(modes))) {
        Extrema <- .findExtrema(modes[, k])
        OUT$count[k] <- abs(nrow(Extrema$maxima) + nrow(Extrema$minima) - Extrema$crossings) <= 1L
        Envelopes <- .buildEnvelopes(modes[, k], extrema = Extrema, boundary = boundary)
        if (is.null(Envelopes)) next
        Deviation <- .measureDeviation(mean = (Envelopes$upper + Envelopes$lower) / 2,
                                       amplitude = (Envelopes$upper - Envelopes$lower) / 2)
        OUT$deviation[k] <- .meetsDoubleThreshold(Deviation, thresholds = thresholds,
                                                  fraction = fraction)
    }
    OUT
}
