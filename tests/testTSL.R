library(xemd)

# Contract of the record operations on a long table of two records with two
# components: decomposition, metrics, selection and reconstruction.
rejects <- function(expr) inherits(tryCatch(expr, error = identity), "error")
isClose <- function(value, reference) {
    all(abs(value - reference) <= 1e-10 * (1 + max(abs(reference))))
}
Time <- seq(0, 9, length.out = 500)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) + 0.5 * Time
TSL <- data.table::rbindlist(lapply(c("AT", "VT"), function(id) {
    data.table::rbindlist(lapply(c("H1", "UP"), function(ocid) {
        data.table::data.table(ID = id, OCID = ocid, t = Time, s = Signal)
    }))
}))
Keys <- c("OCID", "ID")

# listDecompositions names exported methods, the scalar ones first.
Methods <- listDecompositions()
stopifnot(is.character(Methods), length(Methods) == 17L,
          all(vapply(Methods, exists, logical(1L), envir = asNamespace("xemd"))),
          identical(Methods[1L], "decomposeEMD"), "decomposeMVMD" %in% Methods[14:17])

# decomposeTSL with a scalar method: one table per field, keyed by series,
# each series with its own modes, and the identity under a declared grouping.
Decomposition <- decomposeTSL(TSL, method = "decomposeEMD",
                              control = list(sifting = list(stopRule = "doubleThreshold")),
                              by = character())
stopifnot(all(c("signal", "modes", "residue", "diagnostics", "termination") %in% names(Decomposition)),
          all(vapply(Decomposition, data.table::is.data.table, logical(1L))),
          identical(names(Decomposition$modes), c(Keys, "mode", "t", "s")),
          identical(names(Decomposition$signal), c(Keys, "t", "s")),
          nrow(Decomposition$signal) == nrow(TSL),
          nrow(unique(Decomposition$termination[, Keys, with = FALSE])) == 4L,
          all(Decomposition$termination$stopReason %in%
                  c("terminalResidue", "modeLimit", "insufficientExtrema")))
Grouped <- decomposeTSL(data.table::data.table(RecordID = "R1", TSL), method = "decomposeEMD",
                        control = list(sifting = list(stopRule = "doubleThreshold")), by = "RecordID")
stopifnot(identical(names(Grouped$modes), c("RecordID", Keys, "mode", "t", "s")),
          identical(Grouped$modes$s, Decomposition$modes$s))

# An ensemble method draws its noise from the caller's stream by series; a
# coupled method decomposes the components of a record together.
set.seed(628)
Ensemble <- decomposeTSL(TSL, method = "decomposeEEMD",
                         control = list(realizationCount = 2, amplitudeRatio = 0.2,
                                        divisor = "realizations",
                                        sifting = list(stopRule = "doubleThreshold")),
                         by = character())
stopifnot("meanNoise" %in% names(Ensemble), nrow(Ensemble$meanNoise) == nrow(TSL),
          identical(names(Ensemble$termination), c(Keys, "realization", "mode", "iterations", "stopReason")),
          rejects(decomposeTSL(data.table::data.table(realization = 1L, TSL), method = "decomposeEEMD",
                               control = list(realizationCount = 2, amplitudeRatio = 0.2,
                                              divisor = "realizations",
                                              sifting = list(stopRule = "doubleThreshold")),
                               by = "realization")))
Coupled <- decomposeTSL(TSL[TSL$ID == "AT"], method = "decomposeMVMD",
                        control = list(initialCenters = c(0, 1, 3), alpha = 2000 * (9 / 499)^2 / (8 * pi^2),
                                       tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6,
                                       maxIterations = 50),
                        by = character())
stopifnot("fitResidual" %in% names(Coupled), "centers" %in% names(Coupled),
          nrow(Coupled$fitResidual) == 2L * length(Time))
stopifnot(rejects(decomposeTSL(TSL[, c("ID", "t", "s")], method = "decomposeEMD",
                               control = list(sifting = list(stopRule = "doubleThreshold")),
                               by = character())),
          rejects(decomposeTSL(TSL, method = "decomposeNothing", control = list(), by = character())))

# measureModesTSL: metrics by mode and by series, in the order of the tables.
Measures <- measureModesTSL(Decomposition)
stopifnot(identical(names(Measures), c("metrics", "metricsSeries")),
          identical(names(Measures$metrics),
                    c(Keys, "mode", "frequency", "bandwidth", "energy", "energyFraction", "meanPeriod")),
          identical(names(Measures$metricsSeries), c(Keys, "modeCount", "residueFraction", "crossFraction")),
          nrow(Measures$metrics) == nrow(unique(Decomposition$modes[, c(Keys, "mode"), with = FALSE])),
          nrow(Measures$metricsSeries) == 4L, all(Measures$metrics$frequency > 0),
          all(Measures$metrics$energyFraction >= 0), all(Measures$metrics$energyFraction <= 1))

# selectModesTSL and reconstructTSL: keeping every mode with the residue
# reproduces the signal; a band keeps a subset; a listed removal applies to
# the named series only; an incomplete selection is rejected.
Every <- selectModesTSL(Decomposition, metrics = Measures$metrics,
                        criterion = list(rule = "energyFraction", minimum = 0, action = "keep"))
stopifnot(identical(names(Every), c("modes", "selection")),
          identical(names(Every$selection), c(Keys, "mode", "criterion", "retained")),
          all(Every$selection$retained), identical(Every$modes, Decomposition$modes))
Complete <- reconstructTSL(Decomposition, selection = Every$selection)
stopifnot(identical(names(Complete), "reconstruction"),
          identical(names(Complete$reconstruction), c(Keys, "t", "s")),
          isClose(Complete$reconstruction$s, Decomposition$signal$s))
Band <- selectModesTSL(Decomposition, metrics = Measures$metrics,
                       criterion = list(rule = "band", fMin = 2, fMax = 15, action = "keep"))
stopifnot(sum(Band$selection$retained) < nrow(Band$selection), nrow(Band$modes) < nrow(Every$modes),
          identical(unique(Band$selection$criterion), "keep frequency in [2, 15] Hz"))
Filtered <- reconstructTSL(Decomposition, selection = Band$selection, includeResidue = FALSE)
stopifnot(nrow(Filtered$reconstruction) == nrow(TSL),
          !isClose(Filtered$reconstruction$s, Decomposition$signal$s))
Listed <- selectModesTSL(Decomposition, metrics = Measures$metrics,
                         criterion = list(rule = "modes", modes = data.frame(ID = "VT", mode = 1),
                                          action = "remove"))
stopifnot(all(Listed$selection$retained[Listed$selection$ID == "AT"]),
          !all(Listed$selection$retained[Listed$selection$ID == "VT"]))
stopifnot(rejects(reconstructTSL(Decomposition, selection = Every$selection[-1L])),
          rejects(reconstructTSL(Decomposition, selection = Every$selection, includeResidue = NA)),
          rejects(selectModesTSL(Decomposition, metrics = Measures$metrics,
                                 criterion = list(rule = "band", action = "keep"))))
