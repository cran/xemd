# Perturbation realizations w^(m), m = 1, ..., realizationCount, of
# eq-eemd-blancura (eemd.model:27-33; ceemdan.model:29): independent Gaussian
# samples of zero mean and unit variance drawn by stats::rnorm from the
# caller's RNG stream, without a seed of their own. Column m is realization
# m, so the samples are drawn realization by realization.
.generateNoise <- function(sampleCount, realizationCount) {
    isCount <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
    if (!(isCount(sampleCount) && isCount(realizationCount))) {
        stop("sampleCount and realizationCount must be whole numbers of at least one")
    }
    matrix(stats::rnorm(sampleCount * realizationCount), nrow = sampleCount, ncol = realizationCount)
}
