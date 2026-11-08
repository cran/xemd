# Quality of reconstruction factor of an estimate against its reference
# (eq-validacion-qrf), in dB; Inf when both coincide.
measureQRF <- function(estimate, reference) {
    if (!is.numeric(estimate) || !is.null(dim(estimate)) || !all(is.finite(estimate)) ||
        !is.numeric(reference) || !is.null(dim(reference)) ||
        !all(is.finite(reference)) || length(estimate) != length(reference)) {
        stop("'estimate' and 'reference' must be finite numeric vectors of equal length")
    }
    estimate <- as.vector(estimate)
    reference <- as.vector(reference)
    Norm <- sqrt(sum(reference^2))
    if (Norm == 0) stop("'reference' must have a nonzero norm")
    20 * log10(Norm / sqrt(sum((reference - estimate)^2)))
}
