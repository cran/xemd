library(xemd)

# Contract of the variational family: modes, fit residual and centers in Hz,
# the exact fit identity, and the rejection of an invalid configuration.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
SamplingInterval <- 1 / 1000
Time <- seq_len(500) * SamplingInterval
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
Signal <- rowSums(Tones)
Alpha <- 2000 * SamplingInterval^2 / (8 * pi^2)
Fields <- c("modes", "fitResidual", "centers", "iterations", "converged", "stopReason")

OUT <- decomposeVMD(Signal, samplingInterval = SamplingInterval,
                    initialCenters = (0:2) / (6 * SamplingInterval), alpha = Alpha, tau = 0.5,
                    fixedZeroCenter = FALSE, tolerance = 1e-6, maxIterations = 98)
stopifnot(identical(names(OUT), Fields), identical(dim(OUT$modes), c(500L, 3L)),
          length(OUT$fitResidual) == 500L, length(OUT$centers) == 3L,
          all(OUT$centers >= 0), all(OUT$centers <= 1 / (2 * SamplingInterval)),
          identical(OUT$iterations, rep(OUT$iterations[1L], 3L)), length(OUT$converged) == 3L,
          OUT$stopReason %in% c("tolerance", "iterationLimit"),
          isClose(rowSums(OUT$modes) + OUT$fitResidual, Signal),
          rejects(decomposeVMD(Signal, samplingInterval = SamplingInterval, initialCenters = 1:3,
                               alpha = -1, tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6,
                               maxIterations = 98)))

Multichannel <- cbind(Signal, Tones[, 1L], Tones[, 2L])
OUT <- decomposeMVMD(Multichannel, samplingInterval = SamplingInterval,
                     initialCenters = (0:2) / (6 * SamplingInterval), alpha = Alpha, tau = 0.5,
                     fixedZeroCenter = FALSE, tolerance = 1e-6, maxIterations = 98)
stopifnot(identical(names(OUT), Fields), identical(dim(OUT$modes), c(500L, 3L, 3L)),
          identical(dim(OUT$fitResidual), c(500L, 3L)), length(OUT$centers) == 3L,
          OUT$stopReason %in% c("tolerance", "iterationLimit"),
          isClose(apply(OUT$modes, c(1, 3), sum) + OUT$fitResidual, Multichannel),
          rejects(decomposeMVMD(Signal, samplingInterval = SamplingInterval,
                                initialCenters = (0:2) / (6 * SamplingInterval), alpha = Alpha,
                                tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6)))

OUT <- decomposeSVMD(Signal, samplingInterval = SamplingInterval, initialCenter = 0,
                     alpha = 10 * SamplingInterval^2 / (4 * pi^2),
                     maxAlpha = 20000 * SamplingInterval^2 / (4 * pi^2), tau = 0, tolerance = 1e-6,
                     maxIterations = 50, stopping = list(stopRule = "modeCount", modeCount = 3))
stopifnot(identical(names(OUT), Fields), nrow(OUT$modes) == 500L, ncol(OUT$modes) == 3L,
          length(OUT$centers) == 3L, length(OUT$iterations) == 3L, length(OUT$converged) == 3L,
          OUT$stopReason %in% c("modeCount", "noisePower", "modeLimit"),
          isClose(rowSums(OUT$modes) + OUT$fitResidual, Signal),
          rejects(decomposeSVMD(Signal, samplingInterval = SamplingInterval, alpha = 1, maxAlpha = 2,
                                tau = 0, tolerance = 1e-6, stopping = list(stopRule = "never"))))
