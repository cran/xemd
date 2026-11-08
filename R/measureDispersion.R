# Dispersion of the ensemble members around their mean, sample by sample and
# mode by mode (modos.validacion:59-91): the sample variance
# s^2 = sum_m (C^(m) - mean C)^2 / (M - 1) over the M realizations
# (eq-validacion-varianza-muestral) and the variance of the mean s^2 / M.
# realizations is the array samples x modes x realizations of
# .alignRealizations, whose absent modes are zero.
measureDispersion <- function(realizations) {
    if (!is.array(realizations) || length(dim(realizations)) != 3L || !is.numeric(realizations) ||
        !all(is.finite(realizations)) || dim(realizations)[3L] < 2L) {
        stop("'realizations' must be a finite numeric array samples x modes x realizations ",
             "with at least two realizations")
    }
    attributes(realizations) <- list(dim = dim(realizations))
    M <- dim(realizations)[3L]
    Mean <- rowMeans(realizations, dims = 2L)
    Variance <- rowSums((realizations - as.vector(Mean))^2, dims = 2L) / (M - 1)
    list(variance = Variance, meanVariance = Variance / M)
}
