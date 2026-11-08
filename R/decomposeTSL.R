# Record operation on the canonical long table of gmsp, TSL: t in seconds,
# s, ID (the quantity: AT, VT, DT or another) and OCID (the component) with
# metadata (gmsp tableContracts.R:8-11); by names the metadata columns that
# identify a record. A group is a record and an ID; its OCID share one time
# grid. Groups follow ascending by and the IDs AT, VT and DT first, the
# preferred order of gmsp (tableContracts.R:105-108), then the other IDs
# ascending; the OCID of a group ascend; the order of the rows of TSL does
# not intervene. A grid that is not uniform is regularized by
# .regularizeSeries with one warning per group.
# method is one of listDecompositions(); control is the named list of its
# arguments, common to every group, except those .decomposeGroup derives:
# signal, samplingInterval, noise (from realizationCount) and noiseChannels
# (from noiseChannelCount and noiseRatio). AT, VT and DT have different units,
# so the controls in the units of the signal stay in the methods: the rule
# "pair" of the masks and the stops of SVMD other than "modeCount" are
# rejected. Returns data.tables keyed by the by columns, the other metadata,
# OCID and ID: signal, the decomposed series on the output grid (t, s); modes
# (mode, t, s); residue, fitResidual and meanNoise (t, s); diagnostics (mode,
# iterations, converged); termination (the native termination table, or
# stopReason with innerStopReason in itEMD); centers and maskFrequencies
# (mode, frequency in Hz); boundaries
# (boundary, frequency in Hz); each when the outputs of the method carry it.
# Rows follow by, OCID, ID, mode or boundary and t, the order of the tables
# of gmsp (processSignals.R:11-12); every OCID keeps its own number of modes,
# without padding. An output whose samples differ from the grid of its group
# is a defect of the method that stops the operation. TSL is not modified.
decomposeTSL <- function(TSL, method, control, by) {
    DATA <- c("t", "s", "ID", "OCID")
    if (!(is.data.frame(TSL) && nrow(TSL) >= 1L && !anyDuplicated(names(TSL)) &&
          all(DATA %in% names(TSL)))) {
        stop("TSL must be a table with rows, unique column names and the columns t, s, ID and OCID")
    }
    if (!(is.numeric(TSL[["t"]]) && all(is.finite(TSL[["t"]])) &&
          is.numeric(TSL[["s"]]) && all(is.finite(TSL[["s"]])))) {
        stop("t and s must be finite real numbers")
    }
    if (!(is.character(TSL[["ID"]]) && !anyNA(TSL[["ID"]]) &&
          is.character(TSL[["OCID"]]) && !anyNA(TSL[["OCID"]]))) {
        stop("ID and OCID must be character without NA")
    }
    if (!(is.character(by) && !anyNA(by) && !anyDuplicated(by) &&
          all(by %in% setdiff(names(TSL), DATA)))) {
        stop("by must name metadata columns of TSL, or be character()")
    }
    if (!(is.character(method) && length(method) == 1L && method %in% listDecompositions())) {
        stop("method must be one of listDecompositions()")
    }
    META <- setdiff(names(TSL), c(DATA, by))
    COLS <- intersect(c(by, META), c("mode", "iterations", "converged", "stopReason", "innerStopReason",
                                     "frequency", "boundary",
                                     switch(method, decomposeEEMD = "realization",
                                            decomposeCEEMD = c("realization", "sign"),
                                            decomposeMaskedEMD = "phase", decomposeBEMD = "direction")))
    if (length(COLS)) stop("metadata take the names of output columns: ", paste(COLS, collapse = ", "))
    Formals <- names(formals(method))
    Accepted <- c(setdiff(Formals, c("signal", "samplingInterval", "noise", "noiseChannels")),
                  if ("noise" %in% Formals) "realizationCount",
                  if ("noiseChannels" %in% Formals) c("noiseChannelCount", "noiseRatio"))
    if (!(is.list(control) && length(names(control)) == length(control) &&
          all(names(control) %in% Accepted) && !anyDuplicated(names(control)))) {
        stop("control must be a named list of arguments of ", method, " among ",
             paste(Accepted, collapse = ", "))
    }
    Ratio <- control[["noiseRatio"]]
    if ("noiseChannels" %in% Formals &&
        !(is.numeric(Ratio) && length(Ratio) == 1L && is.finite(Ratio) && Ratio >= 0)) {
        stop("control$noiseRatio must be a finite number not below zero")
    }
    if (method == "decomposeMaskedEMD" && is.list(control[["masking"]]) &&
        identical(control[["masking"]][["rule"]], "pair")) {
        stop("masking$rule pair has its amplitude in the units of the signal, which AT, VT and DT do not share")
    }
    if (method == "decomposeSVMD" && is.list(control[["stopping"]]) &&
        !identical(control[["stopping"]][["stopRule"]], "modeCount")) {
        stop("stopping$stopRule must be modeCount: noisePower is in the units of the signal, ",
             "which AT, VT and DT do not share")
    }
    DT <- as.data.table(TSL)
    if (anyDuplicated(DT, by = c(by, "ID", "OCID", "t"))) stop("a time repeats within a series of TSL")
    if (data.table::uniqueN(DT, by = c(by, "ID", "OCID", META)) !=
        data.table::uniqueN(DT, by = c(by, "ID", "OCID"))) {
        stop("the metadata other than by must be constant within each record, ID and OCID")
    }
    AUX <- unique(DT[["ID"]])
    AUX <- c(intersect(c("AT", "VT", "DT"), AUX), sort(setdiff(AUX, c("AT", "VT", "DT")), method = "radix"))
    set(DT, j = "ID", value = factor(DT[["ID"]], levels = AUX))
    data.table::setorderv(DT, c(by, "ID", "OCID", "t"), na.last = TRUE)
    Groups <- lapply(split(DT, by = c(by, "ID"), sorted = FALSE), function(Group) {
        Label <- paste(c(by, "ID"), vapply(c(by, "ID"), function(Column) as.character(Group[[Column]][1L]), ""),
                       sep = " = ", collapse = ", ")
        Times <- split(Group[["t"]], factor(Group[["OCID"]], levels = unique(Group[["OCID"]])))
        if (!all(vapply(Times, identical, logical(1L), Times[[1L]]))) {
            stop(Label, ": the OCID do not share one time grid")
        }
        if (length(Times[[1L]]) < 2L) stop(Label, ": a series needs two samples or more")
        Grid <- .regularizeSeries(time = as.double(Times[[1L]]),
                                  series = matrix(as.double(Group[["s"]]), ncol = length(Times)))
        if (Grid$regularized) {
            warning(warningCondition(paste0(Label, ": the time grid is not uniform; its ", length(Times[[1L]]),
                                            " samples were interpolated by monoH.FC onto ", length(Grid$time),
                                            " samples every ", Grid$samplingInterval, " s"),
                                     class = c("xemdWarning", "moduleWarning"), call = NULL))
        }
        c(Grid, list(keys = unique(Group, by = "OCID")[, c(by, META, "OCID", "ID"), with = FALSE]))
    })
    Parts <- unlist(lapply(Groups, function(Group) {
        Outputs <- .decomposeGroup(series = Group$series, samplingInterval = Group$samplingInterval,
                                   method = method, control = control)
        N <- length(Group$time)
        lapply(seq_along(Outputs), function(c) {
            rows <- function(n) Group$keys[rep(c, n)]
            AUX <- Outputs[[c]]
            K <- ncol(AUX$modes)
            Fields <- intersect(c("residue", "fitResidual", "meanNoise"), names(AUX))
            if (any(c(nrow(AUX$modes), lengths(AUX[Fields])) != N)) {
                stop("defect of ", method, ": an output has a number of samples other than the ", N,
                     " of its grid")
            }
            OUT <- list(signal = data.table(rows(N), t = Group$time, s = Group$series[, c]),
                        modes = data.table(rows(N * K), mode = rep(seq_len(K), each = N),
                                           t = rep(Group$time, K), s = as.vector(AUX$modes)))
            for (Field in Fields) OUT[[Field]] <- data.table(rows(N), t = Group$time, s = AUX[[Field]])
            if (!is.null(AUX$iterations)) {
                OUT$diagnostics <- data.table(rows(K), mode = seq_len(K), iterations = AUX$iterations,
                                              converged = AUX$converged)
            }
            if (!is.null(AUX$termination)) {
                OUT$termination <- data.table(rows(nrow(AUX$termination)), AUX$termination)
            } else if (!is.null(AUX$stopReason)) {
                OUT$termination <- data.table(rows(1L), stopReason = AUX$stopReason,
                                              innerStopReason = AUX$innerStopReason)
            }
            for (Field in intersect(c("centers", "maskFrequencies"), names(AUX))) {
                OUT[[Field]] <- data.table(rows(length(AUX[[Field]])), mode = seq_along(AUX[[Field]]),
                                           frequency = AUX[[Field]])
            }
            if (!is.null(AUX$boundaries)) {
                OUT$boundaries <- data.table(rows(length(AUX$boundaries)), boundary = seq_along(AUX$boundaries),
                                             frequency = AUX$boundaries)
            }
            OUT
        })
    }), recursive = FALSE)
    lapply(stats::setNames(nm = names(Parts[[1L]])), function(Name) {
        AUX <- rbindlist(lapply(Parts, `[[`, Name))
        data.table::setorderv(AUX, intersect(c(by, "OCID", "ID", "mode", "boundary", "t"), names(AUX)),
                              na.last = TRUE)
        set(AUX, j = "ID", value = as.character(AUX[["ID"]]))
        AUX
    })
}
