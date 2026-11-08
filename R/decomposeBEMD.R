# Bivariate empirical mode decomposition (memd.model:9-62). The two columns of
# signal are the real and the imaginary parts of z = x1 + i x2, projected on
# directionCount directions equally spaced on the circle,
# phi_d = 2 pi (d - 1) / D, the grid of Rlibeemd bemd (bemd.implementation:58).
# envelope selects the mean of the sifting: "trajectory", the envelopes of the
# trajectory itself at the maxima of each projection, averaged (scheme 1,
# memd.model:26-36), or "projection", the envelopes of the projections
# oriented on their directions, summed with the factor 2 / D (scheme 2,
# :40-52). Each extraction is the fixed sifting (emd.model:220-226); the
# trajectory envelope admits only reflection. Default: 64 directions, the
# default of Rlibeemd bemd and of emdr memd for the same control. The
# decomposition runs on the bare samples, so no attribute of signal reaches
# the outputs.
decomposeBEMD <- function(signal, directionCount = 64, envelope, sifting) {
    if (!(is.matrix(signal) && is.numeric(signal) && ncol(signal) == 2L && all(is.finite(signal)))) {
        stop("signal must be a numeric matrix of finite real values with two columns")
    }
    if (!(is.numeric(directionCount) && length(directionCount) == 1L && is.finite(directionCount) &&
          directionCount >= 1 && directionCount == round(directionCount))) {
        stop("directionCount must be a whole number of at least one")
    }
    if (!(is.character(envelope) && length(envelope) == 1L && envelope %in% c("trajectory", "projection"))) {
        stop("envelope must be trajectory or projection")
    }
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule != "fixed") stop("sifting$stopRule must be fixed for BEMD")
    if (envelope == "trajectory" && sifting$boundary != "reflection") {
        stop("sifting$boundary must be reflection for the trajectory envelope")
    }
    Angle <- 2 * pi * (seq_len(directionCount) - 1) / directionCount
    .decomposeVector(matrix(as.vector(signal), nrow = nrow(signal)),
                     directions = cbind(cos(Angle), sin(Angle)), envelope = envelope, sifting = sifting,
                     acceptIncomplete = FALSE)
}
