# Mask frequency of the 2005 masking (eq-mascaras-frecuencia-original,
# mascaras.model:52-60): nu = sum(A_1 nu_1^2) / sum(A_1 nu_1), with A_1 and nu_1
# the instantaneous amplitude and frequency of the first IMF without mask,
# y_1 = E_1(signal). The sums run over the samples where nu_1 is defined
# (imf.model:124). Returns a scalar in Hz from samplingInterval in seconds.
estimateMaskFrequency <- function(signal, samplingInterval, sifting) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal))) {
        stop("signal must be a numeric vector of finite real values without dim")
    }
    if (!(is.numeric(samplingInterval) && length(samplingInterval) == 1L &&
          is.finite(samplingInterval) && samplingInterval > 0)) {
        stop("samplingInterval must be a finite number of seconds greater than zero")
    }
    sifting <- .checkSifting(sifting)
    if (!sifting$stopRule %in% c("SD", "Snumber", "doubleThreshold", "fixed")) {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    First <- .extractIMF(as.vector(signal), sifting = sifting)
    if (is.null(First)) {
        stop("signal has no first intrinsic mode function, so its mask frequency is undefined")
    }
    Attributes <- computeInstantaneousAttributes(matrix(First$mode, ncol = 1L),
                                                 samplingInterval = samplingInterval)
    Defined <- !is.na(Attributes$frequency)
    Weight <- Attributes$amplitude[Defined] * Attributes$frequency[Defined]
    if (sum(Weight) == 0) {
        stop(errorCondition("the weights of the mask frequency sum to zero, so it is undefined",
                            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    sum(Weight * Attributes$frequency[Defined]) / sum(Weight)
}
