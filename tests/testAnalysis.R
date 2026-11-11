library(xemd)

# Contract of the mode measures, tests and selections on three known tones.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time), 0.0625 * cos(2 * pi * 288 * Time))
Signal <- rowSums(Tones)

# computeInstantaneousAttributes and buildHilbertSpectrum: the frequency of
# each tone is recovered away from the edges; the spectrum has one column per
# cell and its marginal one value per cell.
Attributes <- computeInstantaneousAttributes(Tones, samplingInterval = SamplingInterval)
Inner <- 101:900
stopifnot(identical(names(Attributes), c("amplitude", "phase", "frequency")),
          identical(dim(Attributes$frequency), dim(Tones)),
          all(abs(Attributes$frequency[Inner, 1L] - 2) < 0.5),
          all(abs(Attributes$frequency[Inner, 2L] - 24) < 1),
          all(abs(Attributes$amplitude[Inner, 2L] - 0.25) < 0.05))
Spectrum <- buildHilbertSpectrum(Attributes$frequency, Attributes$amplitude,
                                 edges = c(0, 10, 100, 500), exponent = 2,
                                 samplingInterval = SamplingInterval)
stopifnot(identical(names(Spectrum), c("spectrum", "marginal")),
          identical(dim(Spectrum$spectrum), c(1000L, 3L)), length(Spectrum$marginal) == 3L,
          which.max(Spectrum$marginal) == 1L,
          rejects(buildHilbertSpectrum(Attributes$frequency, Attributes$amplitude,
                                       edges = c(10, 0), exponent = 2,
                                       samplingInterval = SamplingInterval)))

# Measures: orthogonality of the tones, the index between adjacent modes,
# the quality in dB, the energy order, and one factor or distance per mode.
Orthogonality <- measureOrthogonality(Tones)
stopifnot(identical(names(Orthogonality), c("innerProduct", "cosine")),
          identical(dim(Orthogonality$cosine), c(3L, 3L)),
          all(abs(Orthogonality$cosine[upper.tri(Orthogonality$cosine)]) < 0.01),
          all(abs(diag(Orthogonality$cosine) - 1) < 1e-12))
stopifnot(length(measurePMSI(Tones, denominator = "homogeneous")) == 2L,
          length(measurePMSI(Tones, denominator = "regularized")) == 2L,
          rejects(measurePMSI(Tones, denominator = "other")))
stopifnot(measureQRF(estimate = Signal, reference = Signal) >
              measureQRF(estimate = Signal - Tones[, 3L], reference = Signal),
          measureQRF(estimate = Signal - Tones[, 3L], reference = Signal) > 20)
stopifnot(identical(selectModesEnergy(Tones), 1L),
          identical(selectModesEnergy(Tones[, c(3, 1, 2)]), 2L),
          rejects(selectModesEnergy(cbind(Tones[, 1L], Tones[, 1L]))))
stopifnot(length(measureFCF(signal = Signal, modes = Tones)) == 3L,
          length(measureDTW(signal = Signal, modes = Tones)) == 3L,
          all(measureDTW(signal = Signal, modes = Tones) >= 0),
          rejects(measureFCF(signal = Signal[-1L], modes = Tones)))

# Tests of the modes: the IMF conditions, the white-noise bounds and the
# fractional-noise thresholds, one value per column.
IMF <- assessIMF(Tones)
stopifnot(identical(names(IMF), c("count", "deviation")), length(IMF$count) == 3L,
          is.logical(IMF$count), length(IMF$deviation) == 3L, all(IMF$count))
Scaled <- Tones / sqrt(mean((Signal - mean(Signal))^2))
for (Bound in list(list(bound = "gaussian", rescalingMode = NULL),
                   list(bound = "chiSquare", rescalingMode = 1))) {
    White <- assessWhiteNoise(Scaled, probabilities = c(0.025, 0.975), bound = Bound$bound,
                              rescalingMode = Bound$rescalingMode)
    stopifnot(identical(names(White), c("energy", "period", "lower", "upper", "outside")),
              length(White$outside) == 3L, is.logical(White$outside), all(White$period > 0),
              all(White$lower <= White$upper))
}
stopifnot(rejects(assessWhiteNoise(Scaled, probabilities = c(0.975, 0.025), bound = "gaussian",
                                   rescalingMode = NULL)))
Long <- seq_len(2048) / 1000
Fractional <- assessFractionalNoise(cbind(cos(2 * pi * 2 * Long), 0.25 * cos(2 * pi * 24 * Long),
                                          0.0625 * cos(2 * pi * 288 * Long)),
                                    hurst = 0.5, probability = 0.95)
stopifnot(identical(names(Fractional), c("energy", "reference", "upper", "outside")),
          length(Fractional$energy) == 3L, is.na(Fractional$outside[1L]),
          is.logical(Fractional$outside),
          rejects(assessFractionalNoise(Tones, hurst = 0.5, probability = 0.95)))

# measureDispersion: the variance by sample and mode across realizations.
Realizations <- array(c(2, 4, 6, 8, 1, 1, 1, 1, 4, 6, 8, 10, 0, 0, 0, 0), dim = c(4, 2, 2))
Dispersion <- measureDispersion(Realizations)
stopifnot(identical(names(Dispersion), c("variance", "meanVariance")),
          identical(dim(Dispersion$variance), c(4L, 2L)),
          isTRUE(all.equal(Dispersion$variance[, 1L], rep(2, 4L))),
          isTRUE(all.equal(Dispersion$meanVariance, Dispersion$variance / 2)),
          rejects(measureDispersion(Tones)))
