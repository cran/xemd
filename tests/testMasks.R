library(xemd)

# Contract of the masked decompositions and of TVF-EMD on a toy case.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
Time <- seq(0, 9, length.out = 500)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) + 0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", maxModes = 3)
Fields <- c("modes", "residue", "iterations", "converged", "stopReason")

# estimateMaskFrequency: a frequency in the units of the sampling interval.
Frequency <- estimateMaskFrequency(Signal, samplingInterval = 1, sifting = Sifting)
stopifnot(is.numeric(Frequency), length(Frequency) == 1L, Frequency > 0, Frequency < 0.5,
          isClose(estimateMaskFrequency(Signal, samplingInterval = 9 / 499, sifting = Sifting),
                  Frequency * 499 / 9))

# decomposeMaskedEMD: the uniform-phase rule and the pair rule, each with one
# mask frequency per mode; an unknown rule is rejected.
OUT <- decomposeMaskedEMD(Signal, samplingInterval = 1,
                          masking = list(rule = "uniformPhase", amplitudeRatio = 1, phaseCount = 4),
                          sifting = Sifting)
stopifnot(identical(names(OUT), c(Fields, "maskFrequencies", "termination")),
          nrow(OUT$modes) == length(Signal), ncol(OUT$modes) <= 3L,
          length(OUT$maskFrequencies) == ncol(OUT$modes),
          OUT$stopReason %in% c("terminalResidue", "modeLimit", "insufficientExtrema", "incompleteExtraction"),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal))
stopifnot(identical(names(OUT$termination), c("mode", "phase", "iterations", "stopReason")),
          nrow(OUT$termination) >= 1L)
Pair <- decomposeMaskedEMD(Signal, samplingInterval = 1,
                           masking = list(rule = "pair", frequency = Frequency,
                                          amplitude = sqrt(mean((Signal - mean(Signal))^2))),
                           sifting = Sifting)
stopifnot(identical(names(Pair), names(OUT)), isClose(rowSums(Pair$modes) + Pair$residue, Signal),
          rejects(decomposeMaskedEMD(Signal, samplingInterval = 1, masking = list(rule = "random"),
                                     sifting = Sifting)))

# decomposeIteratedMaskedEMD: the last masked decomposition with the
# frequencies that produced it and the outer termination.
OUT <- decomposeIteratedMaskedEMD(Signal, samplingInterval = 1,
                                  initialFrequencies = 55 / 500 / 2 / 2^(0:2), phaseCount = 4,
                                  amplitudeRatio = 1, tolerance = 0.1, maxIterations = 15,
                                  sifting = Sifting)
stopifnot(identical(names(OUT), c("modes", "residue", "iterations", "converged", "stopReason",
                                  "maskFrequencies", "maskIterations", "innerStopReason")),
          nrow(OUT$modes) == length(Signal), OUT$stopReason %in% c("tolerance", "iterationLimit"),
          length(OUT$maskFrequencies) == ncol(OUT$modes), is.integer(OUT$maskIterations),
          OUT$maskIterations >= 1L, OUT$maskIterations <= 15L,
          OUT$innerStopReason %in% c("terminalResidue", "modeLimit", "insufficientExtrema"),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          rejects(decomposeIteratedMaskedEMD(Signal, samplingInterval = 1, initialFrequencies = -1,
                                             phaseCount = 4, amplitudeRatio = 1, sifting = Sifting)))

# decomposeTVFEMD: components and residue reconstruct the signal.
OUT <- decomposeTVFEMD(Signal, bandwidthThreshold = 0.1, splineOrder = 28,
                       intermittencyThreshold = 0.25, maxSift = 5, maxModes = 3)
stopifnot(identical(names(OUT), Fields), nrow(OUT$modes) == length(Signal), ncol(OUT$modes) <= 3L,
          length(OUT$iterations) == ncol(OUT$modes), OUT$stopReason %in% c("terminalResidue", "modeLimit"),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          rejects(decomposeTVFEMD(Signal, bandwidthThreshold = -0.1, splineOrder = 28,
                                  intermittencyThreshold = 0.25, maxSift = 5)))
