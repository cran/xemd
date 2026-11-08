library(xemd)

# Contract of the vector decompositions: arrays by channel, reconstruction
# per channel, the default direction count, and the rejection of a vector
# that is not a matrix of two or more channels.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
Time <- 2 * pi * (0:255) / 256
z <- cos(0.3 * Time) * exp(2i * Time) + 0.3 * abs(sin(2.3 * Time)) * exp(17i * Time)
Signal <- cbind(Re(z), Im(z))
Sifting <- list(stopRule = "fixed", fixedSift = 5, boundary = "reflection", maxModes = 3)
Fields <- c("modes", "residue", "iterations", "converged", "stopReason")
checkVector <- function(OUT, signal, fields) {
    stopifnot(identical(names(OUT), fields), length(dim(OUT$modes)) == 3L,
              dim(OUT$modes)[1L] == nrow(signal), dim(OUT$modes)[3L] == ncol(signal),
              dim(OUT$modes)[2L] <= 3L, identical(dim(OUT$residue), dim(signal)),
              length(OUT$iterations) == dim(OUT$modes)[2L], all(is.na(OUT$converged)),
              OUT$stopReason %in% c("terminalResidue", "modeLimit", "insufficientExtrema"),
              isClose(apply(OUT$modes, c(1, 3), sum) + OUT$residue, signal))
}

OUT <- decomposeBEMD(Signal, directionCount = 32, envelope = "projection", sifting = Sifting)
checkVector(OUT, Signal, fields = c(Fields, "termination"))
Trajectory <- decomposeBEMD(Signal, directionCount = 32, envelope = "trajectory",
                            sifting = list(stopRule = "fixed", fixedSift = 5, maxModes = 3))
checkVector(Trajectory, Signal, fields = c(Fields, "termination"))
stopifnot(identical(names(OUT$termination), c("mode", "direction", "iterations", "stopReason")),
          nrow(OUT$termination) == 1L, identical(names(Trajectory$termination), names(OUT$termination)))
stopifnot(identical(decomposeBEMD(Signal, envelope = "projection", sifting = Sifting),
                    decomposeBEMD(Signal, directionCount = 64, envelope = "projection",
                                  sifting = Sifting)),
          rejects(decomposeBEMD(Signal, envelope = "complex", sifting = Sifting)),
          rejects(decomposeBEMD(Signal[, 1L], envelope = "projection", sifting = Sifting)))

OUT <- decomposeMEMD(Signal, directionCount = 32, sifting = Sifting)
checkVector(OUT, Signal, fields = Fields)
Three <- decomposeMEMD(cbind(Signal, Signal[, 1L] * Signal[, 2L]), directionCount = 32,
                       sifting = Sifting)
checkVector(Three, cbind(Signal, Signal[, 1L] * Signal[, 2L]), fields = Fields)
stopifnot(rejects(decomposeMEMD(Signal[, 1L], sifting = Sifting)),
          rejects(decomposeMEMD(Signal, directionCount = 0, sifting = Sifting)))

set.seed(628)
Noise <- 0.02 * sqrt(mean(apply(Signal, 2, stats::sd)^2)) * matrix(stats::rnorm(256 * 2), nrow = 256)
RNG <- .Random.seed
OUT <- decomposeNAMEMD(Signal, noiseChannels = Noise, directionCount = 32, sifting = Sifting)
checkVector(OUT, Signal, fields = Fields)
stopifnot(identical(.Random.seed, RNG),
          rejects(decomposeNAMEMD(Signal, noiseChannels = Noise[-1L, ], sifting = Sifting)))
