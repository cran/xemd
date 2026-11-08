library(xemd)

# Contract of the EMD family on a toy case: structure, reconstruction, the
# caller's random stream left untouched, and the rejection of invalid input.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
Time <- seq(0, 9, length.out = 500)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) + 0.5 * Time
Sifting <- list(stopRule = "doubleThreshold")
set.seed(628)
Noise <- matrix(stats::rnorm(3L * length(Signal)), ncol = 3L)
RNG <- .Random.seed

# decomposeEMD: modes and residue reconstruct the signal; the defaults of the
# rule equal the explicit configuration; fixed sifting reports NA convergence.
OUT <- decomposeEMD(Signal, sifting = Sifting)
stopifnot(identical(names(OUT), c("modes", "residue", "iterations", "converged", "stopReason")),
          is.matrix(OUT$modes), nrow(OUT$modes) == length(Signal), ncol(OUT$modes) >= 1L,
          is.null(dimnames(OUT$modes)), length(OUT$residue) == length(Signal),
          is.integer(OUT$iterations), length(OUT$iterations) == ncol(OUT$modes),
          is.logical(OUT$converged), length(OUT$converged) == ncol(OUT$modes),
          OUT$stopReason %in% c("terminalResidue", "modeLimit", "insufficientExtrema"),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          identical(decomposeEMD(Signal, sifting = list(stopRule = "doubleThreshold",
                                                        thresholds = c(0.05, 0.5), fraction = 0.05,
                                                        maxSift = 2000, boundary = "reflection",
                                                        maxModes = Inf)), OUT))
Fixed <- decomposeEMD(Signal, sifting = list(stopRule = "fixed", fixedSift = 10, maxModes = 2))
stopifnot(ncol(Fixed$modes) == 2L, identical(Fixed$iterations, c(10L, 10L)),
          all(is.na(Fixed$converged)), identical(Fixed$stopReason, "modeLimit"))
stopifnot(rejects(decomposeEMD(Signal, sifting = list(stopRule = "consecutive"))),
          rejects(decomposeEMD(cbind(Signal, Signal), sifting = Sifting)))

# decomposeEEMD: the mean of the realizations, with its mean noise; both
# divisors are accepted; the noise must match the signal length.
OUT <- decomposeEEMD(Signal, noise = Noise, amplitudeRatio = 0.2, divisor = "realizations",
                     sifting = Sifting)
stopifnot(identical(names(OUT), c("modes", "residue", "meanNoise", "iterations", "converged", "termination")),
          nrow(OUT$modes) == length(Signal), length(OUT$meanNoise) == length(Signal),
          length(OUT$iterations) == ncol(OUT$modes), length(OUT$converged) == ncol(OUT$modes),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal + OUT$meanNoise),
          identical(.Random.seed, RNG))
stopifnot(identical(names(OUT$termination), c("realization", "mode", "iterations", "stopReason")),
          nrow(OUT$termination) == ncol(Noise))
Available <- decomposeEEMD(Signal, noise = Noise, amplitudeRatio = 0.2, divisor = "available",
                           sifting = Sifting)
stopifnot(identical(dim(Available$modes), dim(OUT$modes)),
          rejects(decomposeEEMD(Signal, noise = Noise[-1L, ], divisor = "realizations",
                                sifting = Sifting)),
          rejects(decomposeEEMD(Signal, noise = Noise, divisor = "mean", sifting = Sifting)))

# decomposeCEEMD: pairs of opposite perturbations, no mean noise.
OUT <- decomposeCEEMD(Signal, noise = Noise, amplitudeRatio = 0.2, sifting = Sifting)
stopifnot(identical(names(OUT), c("modes", "residue", "iterations", "converged", "termination")),
          nrow(OUT$modes) == length(Signal), length(OUT$iterations) == ncol(OUT$modes),
          isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          identical(.Random.seed, RNG),
          rejects(decomposeCEEMD(Signal, noise = Noise[-1L, ], sifting = Sifting)))
stopifnot(identical(names(OUT$termination), c("realization", "sign", "mode", "iterations", "stopReason")),
          nrow(OUT$termination) == 2L * ncol(Noise))

# decomposeCEEMDAN and decomposeICEEMDAN: a common residue by stage, the mode
# limit, and the zero-energy realization signalled with its own class.
Limit <- list(stopRule = "doubleThreshold", maxModes = 3)
for (decompose in list(function(...) decomposeCEEMDAN(...),
                       function(...) decomposeICEEMDAN(..., normalizeNoise = TRUE))) {
    OUT <- decompose(Signal, noise = Noise, amplitudeRatio = 0.2, sifting = Sifting)
    stopifnot(identical(names(OUT), c("modes", "residue", "iterations", "converged", "stopReason")),
              nrow(OUT$modes) == length(Signal), ncol(OUT$modes) >= 2L,
              length(OUT$iterations) == ncol(OUT$modes), length(OUT$converged) == ncol(OUT$modes),
              OUT$stopReason %in% c("terminalResidue", "modeLimit", "insufficientExtrema"),
              isClose(rowSums(OUT$modes) + OUT$residue, Signal),
              identical(.Random.seed, RNG))
    Limited <- decompose(Signal, noise = Noise, amplitudeRatio = 0.2, sifting = Limit)
    stopifnot(ncol(Limited$modes) <= 3L, isClose(rowSums(Limited$modes) + Limited$residue, Signal))
}
# A realization of zero deviation leaves the first scale of CEEMDAN undefined
# and is signalled with its own class; ICEEMDAN, which scales the noise modes,
# gives that realization a null perturbation and keeps it in the divisor.
Silent <- Noise
Silent[, 1L] <- 0
stopifnot(inherits(tryCatch(decomposeCEEMDAN(Signal, noise = Silent, sifting = Sifting),
                            xemdZeroEnergy = identity), "xemdZeroEnergy"))
OUT <- decomposeICEEMDAN(Signal, noise = Silent, normalizeNoise = TRUE, sifting = Sifting)
stopifnot(isClose(rowSums(OUT$modes) + OUT$residue, Signal),
          rejects(decomposeICEEMDAN(Signal, noise = Noise, normalizeNoise = NA, sifting = Sifting)))

# resiftModes: the chain of the modes of an ensemble reconstructs their sum.
Modes <- decomposeEEMD(Signal, noise = Noise, amplitudeRatio = 0.2, divisor = "realizations",
                       sifting = Sifting)$modes
OUT <- resiftModes(Modes, sifting = Sifting)
stopifnot(identical(names(OUT), c("modes", "residue", "iterations", "converged", "stopReason")),
          nrow(OUT$modes) == length(Signal), OUT$stopReason %in% c("chainComplete", "insufficientExtrema"),
          isClose(rowSums(OUT$modes) + OUT$residue, rowSums(Modes)),
          rejects(resiftModes(Signal, sifting = Sifting)))
