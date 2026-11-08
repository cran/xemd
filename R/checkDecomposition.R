# Input of the analysis record operations: decomposition, the list of
# data.tables of decomposeTSL with signal, modes and at most one remainder,
# residue or fitResidual (EWT returns neither). The keys are the columns of
# signal other than t and s, the by columns, the other metadata, OCID and ID;
# they cannot take mode or the output columns of the calling operation,
# columns. Each series of signal has two samples or more at distinct times,
# and every mode and every remainder of a series holds exactly the times of
# that series, so the operations do not depend on the order of the rows.
# Returns the keys.
.checkDecomposition <- function(decomposition, columns) {
    Remainder <- intersect(c("residue", "fitResidual"), names(decomposition))
    if (!(is.list(decomposition) && data.table::is.data.table(decomposition$signal) &&
          all(vapply(decomposition[c("modes", Remainder)], data.table::is.data.table, logical(1L))))) {
        stop("decomposition must be the list of tables of decomposeTSL, with signal and modes")
    }
    if (length(Remainder) > 1L) stop("decomposition must hold residue or fitResidual, not both")
    Tables <- decomposition[c("modes", Remainder)]
    Signal <- decomposition$signal
    Keys <- setdiff(names(Signal), c("t", "s"))
    holds <- function(table, columns) length(table) == length(columns) && setequal(names(table), columns)
    if (!(all(c("t", "s", "OCID", "ID") %in% names(Signal)) && holds(Signal, c(Keys, "t", "s")) &&
          holds(decomposition$modes, c(Keys, "mode", "t", "s")) &&
          all(vapply(decomposition[Remainder], holds, logical(1L), columns = c(Keys, "t", "s"))))) {
        stop("the tables of decomposition must hold the keys of signal, OCID and ID among them, with t and s, ",
             "and mode in modes")
    }
    Taken <- intersect(Keys, c("mode", columns))
    if (length(Taken)) stop("metadata take the names of output columns: ", paste(Taken, collapse = ", "))
    if (!all(vapply(c(list(Signal), Tables), function(Table) {
        is.numeric(Table$t) && all(is.finite(Table$t)) && is.numeric(Table$s) && all(is.finite(Table$s))
    }, logical(1L)))) {
        stop("t and s must be finite real numbers")
    }
    Index <- decomposition$modes$mode
    if (!(is.numeric(Index) && all(is.finite(Index) & Index >= 1 & Index == round(Index)))) {
        stop("mode must hold whole numbers from 1")
    }
    Series <- unique(Signal[, Keys, with = FALSE])
    Samples <- tabulate(Series[Signal, on = Keys, which = TRUE], nbins = nrow(Series))
    if (anyDuplicated(Signal, by = c(Keys, "t")) || any(Samples < 2L)) {
        stop("each series of signal needs two samples or more at distinct times")
    }
    Modes <- unique(decomposition$modes[, c(Keys, "mode"), with = FALSE])
    Lengths <- tabulate(Modes[decomposition$modes, on = c(Keys, "mode"), which = TRUE], nbins = nrow(Modes))
    Complete <- identical(Lengths, Samples[Series[Modes, on = Keys, which = TRUE]]) &&
        all(vapply(decomposition[Remainder], function(Table) {
            identical(tabulate(Series[Table, on = Keys, which = TRUE], nbins = nrow(Series)), Samples)
        }, logical(1L)))
    Outside <- vapply(Tables, function(Table) {
        anyDuplicated(Table, by = c(intersect(c(Keys, "mode"), names(Table)), "t")) > 0L ||
            nrow(Table[!Signal, on = c(Keys, "t")]) > 0L
    }, logical(1L))
    if (!Complete || any(Outside)) stop("the modes and the remainder of each series must lie on its grid in signal")
    Keys
}
