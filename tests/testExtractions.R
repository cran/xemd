library(xemd)

# Exact internal inputs from the six RSN245 cases, without changed sifting,
# phases or directions.
Fixtures <- readRDS("fixtures/extractionFailures.rds")
isClose <- function(value, reference) {
    all(abs(value - reference) <= Fixtures$tolerance * (1 + max(abs(reference))))
}
for (Case in Fixtures$scalar) {
    OUT <- xemd:::.decomposeScalar(Case$signal, sifting = Case$sifting, acceptIncomplete = FALSE)
    stopifnot(ncol(OUT$modes) == 0L, identical(OUT$residue, Case$signal),
              length(OUT$converged) == 0L, OUT$termination$mode == 1L,
              identical(OUT$termination$iterations, Case$iterations),
              OUT$termination$stopReason == "insufficientExtrema")
    # The unchanged EMD entry still exposes its former incomplete candidate.
    Legacy <- decomposeEMD(Case$signal, sifting = Case$sifting)
    stopifnot(isFALSE(Legacy$converged[1L]), identical(Legacy$modes[, 1L], Case$mode),
              isClose(rowSums(Legacy$modes) + Legacy$residue, Case$signal))
}
Case <- Fixtures$initial
OUT <- xemd:::.decomposeScalar(Case$signal, sifting = Case$sifting, acceptIncomplete = FALSE)
stopifnot(ncol(OUT$modes) == 0L, identical(OUT$residue, Case$signal),
          OUT$termination$iterations == 0L, OUT$termination$stopReason == "insufficientExtrema")

for (Case in Fixtures$masked) {
    OUT <- xemd:::.decomposeMasked(Case$signal, schedule = function(k, residue) Case$stage,
                                   phases = Case$phases, sifting = Case$sifting, acceptIncomplete = FALSE)
    stopifnot(ncol(OUT$modes) == 0L, identical(OUT$residue, Case$signal),
              OUT$stopReason == "incompleteExtraction", all(OUT$termination$mode == 1L),
              identical(OUT$termination$phase, Case$phase), all(OUT$termination$iterations == 0L))
}
stopifnot(identical(OUT$termination$stopReason, c("insufficientExtrema", "terminalResidue")))

Case <- Fixtures$vector
OUT <- xemd:::.decomposeVector(Case$signal, directions = Case$directions, envelope = Case$envelope,
                               sifting = Case$sifting, acceptIncomplete = FALSE)
stopifnot(dim(OUT$modes)[2L] == 0L, identical(OUT$residue, Case$signal),
          OUT$termination$mode == 1L, OUT$termination$iterations == 1L,
          OUT$termination$direction == 12L, OUT$termination$stopReason == "insufficientExtrema")
# The failed candidate itself has no envelope at its initial state.
OUT <- xemd:::.decomposeVector(Case$mode, directions = Case$directions, envelope = Case$envelope,
                               sifting = Case$sifting, acceptIncomplete = FALSE)
stopifnot(dim(OUT$modes)[2L] == 0L, identical(OUT$residue, Case$mode),
          OUT$termination$iterations == 0L, OUT$termination$direction == 12L,
          OUT$termination$stopReason == "insufficientExtrema")

# A candidate that reaches the iteration limit is kept, flagged FALSE, and the
# extraction continues: the strict policy stops only on missing envelopes.
Time <- seq(0, 9, length.out = 500)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) + 0.5 * Time
Limit <- list(stopRule = "SD", tolerance = 1e-12, maxSift = 1, boundary = "reflection")
Strict <- xemd:::.decomposeScalar(Signal, sifting = xemd:::.checkSifting(Limit), acceptIncomplete = FALSE)
Legacy <- xemd:::.decomposeScalar(Signal, sifting = xemd:::.checkSifting(Limit), acceptIncomplete = TRUE)
stopifnot(ncol(Strict$modes) >= 1L, all(Strict$converged %in% FALSE), all(Strict$iterations == 1L),
          Strict$termination$stopReason != "iterationLimit", identical(Strict$modes, Legacy$modes))
set.seed(628)
Noise <- matrix(rnorm(500 * 3), 500, 3)
OUT <- decomposeEEMD(Signal, noise = Noise, amplitudeRatio = 0.2, divisor = "realizations", sifting = Limit)
stopifnot(ncol(OUT$modes) >= 1L, all(OUT$converged %in% FALSE),
          all(OUT$termination$stopReason != "iterationLimit"))
OUT <- decomposeMaskedEMD(Signal, samplingInterval = 1,
                          masking = list(rule = "uniformPhase", amplitudeRatio = 1, phaseCount = 4),
                          sifting = Limit)
stopifnot(ncol(OUT$modes) >= 1L, all(OUT$converged %in% FALSE))
