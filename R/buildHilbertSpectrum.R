# Hilbert spectrum H_p as mass per cell (eq-hilbert-espectro and
# eq-hilbert-banda, imf.model:271-289): spectrum[n, b] sums a_k[n]^p over the
# modes k whose frequency at sample n lies in [edges[b], edges[b + 1]), with
# p = exponent, in units of the signal to the power p. A frequency outside the
# edges or NA, where the phase is undefined, adds nothing (imf.model:285). The
# marginal spectrum of each cell is the time integral samplingInterval *
# sum_n spectrum[n, b] (eq-hilbert-marginal-banda, imf.model:299-305), in
# units of the signal to the power p times seconds. No scaling or division by
# the cell width is applied: mass, not density (imf.model:289).
buildHilbertSpectrum <- function(frequency, amplitude, edges, exponent, samplingInterval) {
    if (!is.matrix(frequency) || !is.numeric(frequency) ||
        !all(is.finite(frequency) | is.na(frequency))) {
        stop("'frequency' must be a numeric matrix of finite values or NA")
    }
    if (!is.matrix(amplitude) || !is.numeric(amplitude) ||
        !identical(dim(amplitude), dim(frequency)) || !all(is.finite(amplitude)) ||
        any(amplitude < 0)) {
        stop("'amplitude' must be a finite non-negative matrix with the shape of 'frequency'")
    }
    if (!is.numeric(edges) || !is.null(dim(edges)) || length(edges) < 2L ||
        !all(is.finite(edges)) || any(diff(edges) <= 0)) {
        stop("'edges' must be at least two finite, strictly increasing numbers")
    }
    if (!is.numeric(exponent) || length(exponent) != 1L || !(exponent %in% c(1, 2))) {
        stop("'exponent' must be 1 or 2")
    }
    if (!is.numeric(samplingInterval) || length(samplingInterval) != 1L ||
        !is.finite(samplingInterval) || samplingInterval <= 0) {
        stop("'samplingInterval' must be a finite positive number of seconds")
    }
    attributes(frequency) <- list(dim = dim(frequency))
    attributes(amplitude) <- list(dim = dim(amplitude))
    Spectrum <- matrix(0, nrow = nrow(frequency), ncol = length(edges) - 1L)
    # Within a mode each sample lies in at most one cell, so the masses of a
    # mode never collide and are added to those of the previous modes.
    for (k in seq_len(ncol(frequency))) {
        Cell <- findInterval(frequency[, k], vec = edges)
        IDX <- which(Cell >= 1L & Cell <= ncol(Spectrum))
        AUX <- cbind(IDX, Cell[IDX])
        Spectrum[AUX] <- Spectrum[AUX] + amplitude[IDX, k]^exponent
    }
    list(spectrum = Spectrum, marginal = samplingInterval * colSums(Spectrum))
}
