computeInstantaneousAttributes <- function(modes, samplingInterval) {
    if (!is.matrix(modes) || !is.double(modes) || nrow(modes) < 3L ||
        !all(is.finite(modes))) {
        stop("'modes' must be a finite double matrix with at least 3 rows")
    }
    if (length(samplingInterval) != 1L || !is.numeric(samplingInterval) ||
        !is.finite(samplingInterval) || samplingInterval <= 0) {
        stop("'samplingInterval' must be a finite positive number of seconds")
    }
    attributes(modes) <- list(dim = dim(modes))
    n <- nrow(modes)
    Analytic <- .buildAnalytic(modes)
    Phase <- .unwrapPhase(Analytic)
    # f_c = omega_c / (2 pi) in Hz, with omega_c = phi' (eq-hilbert-frecuencia),
    # by the central difference (phi[n + 1] - phi[n - 1]) / (4 pi dt); undefined
    # where any of the three phases is, and each border repeats its neighbour.
    Frequency <- (Phase[-c(1L, 2L), , drop = FALSE] -
                      Phase[-c(n - 1L, n), , drop = FALSE]) / (4 * pi * samplingInterval)
    Frequency[is.na(Phase[-c(1L, n), , drop = FALSE])] <- NA_real_
    list(amplitude = Mod(Analytic), phase = Phase,
         frequency = Frequency[c(1L, seq_len(n - 2L), n - 2L), , drop = FALSE])
}
