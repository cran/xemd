library(xemd)

# Contract of the empirical wavelet transform, local mean decomposition,
# iterative filtering and the synchrosqueezed transform on three tones.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
SamplingInterval <- 1 / 1000
Time <- seq_len(500) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) + 0.0625 * cos(2 * pi * 288 * Time)
Fields <- c("modes", "residue", "iterations", "converged", "stopReason")

# decomposeEWT: one band per boundary interval, boundaries in Hz below the
# Nyquist frequency, automatic count with a detection list.
OUT <- decomposeEWT(Signal, samplingInterval = SamplingInterval, bandCount = 3)
stopifnot(identical(names(OUT), c("modes", "boundaries")), identical(dim(OUT$modes), c(500L, 3L)),
          length(OUT$boundaries) == 2L, all(diff(OUT$boundaries) > 0),
          all(OUT$boundaries < 1 / (2 * SamplingInterval)),
          isClose(rowSums(OUT$modes), Signal))
Automatic <- decomposeEWT(Signal, samplingInterval = SamplingInterval, bandCount = NULL,
                          detection = list(detector = "locmax", regularization = "none",
                                           thresholdRatio = 0.3))
stopifnot(ncol(Automatic$modes) == length(Automatic$boundaries) + 1L,
          rejects(decomposeEWT(Signal, samplingInterval = SamplingInterval, bandCount = 1)),
          rejects(decomposeEWT(Signal, samplingInterval = SamplingInterval, bandCount = NULL,
                               detection = list(detector = "unknown"))))

# decomposeLMD: product functions and residue reconstruct the signal.
OUT <- decomposeLMD(Signal, maxIterations = 50, envelopeTolerance = 0.01, changeTolerance = 0.01,
                    maxSmoothPasses = 12, maxModes = 3)
stopifnot(identical(names(OUT), Fields), nrow(OUT$modes) == 500L, ncol(OUT$modes) <= 3L,
          length(OUT$iterations) == ncol(OUT$modes), OUT$stopReason %in% c("terminalResidue", "modeLimit"),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          rejects(decomposeLMD(Signal, maxIterations = 0, envelopeTolerance = 0.01,
                               changeTolerance = 0.01, maxSmoothPasses = 12)))

# decomposeIF: the components of a given filter under two boundary rules.
Filter <- c(6.6306270449357867e-06, 0.00053827687799438354, 0.0082858777882808092,
            0.044245208647455882, 0.11927281750086954, 0.20513332487239222, 0.24503568511284984,
            0.20513332487239222, 0.11927281750086952, 0.044245208647455903, 0.0082858777882808109,
            0.00053827687799438213, 6.6306270449358087e-06)
for (Boundary in c("periodic", "zero")) {
    OUT <- decomposeIF(Signal, filter = Filter, boundary = Boundary, tolerance = sqrt(0.001),
                       maxIterations = 50, stretch = 1.6, maxModes = 2)
    stopifnot(identical(names(OUT), Fields), nrow(OUT$modes) == 500L, ncol(OUT$modes) <= 2L,
              length(OUT$iterations) == ncol(OUT$modes),
              OUT$stopReason %in% c("terminalResidue", "modeLimit"),
              isClose(rowSums(OUT$modes) + OUT$residue, Signal))
}
stopifnot(rejects(decomposeIF(Signal, filter = Filter, boundary = "mirror", tolerance = sqrt(0.001),
                              maxIterations = 50, stretch = 1.6)))

# transformSST and extractSSTModes: a complex transform by frequency cell and
# one component per ridge, with the energy of the first tone dominant.
SST <- transformSST(Signal, samplingInterval = SamplingInterval, voiceCount = 8,
                    morletCenter = 2 * pi, threshold = sqrt(.Machine$double.eps))
stopifnot(identical(names(SST), c("transform", "frequency")), is.complex(SST$transform),
          nrow(SST$transform) == 500L, ncol(SST$transform) == length(SST$frequency),
          all(SST$frequency > 0), all(diff(SST$frequency) > 0))
Modes <- extractSSTModes(SST$transform, frequency = SST$frequency, modeCount = 2,
                         reconstructionWidth = 0, morletCenter = 2 * pi)
stopifnot(is.matrix(Modes), identical(dim(Modes), c(500L, 2L)), is.numeric(Modes),
          mean(Modes[, 1L]^2) > mean(Modes[, 2L]^2),
          rejects(transformSST(Signal, samplingInterval = 0, threshold = 0)),
          rejects(extractSSTModes(SST$transform, frequency = SST$frequency, modeCount = 0,
                                  reconstructionWidth = 0, morletCenter = 2 * pi)))
