# Successive extraction of vector modes (memd.model:138-156): r_0 is the
# signal and r_k = r_{k-1} - c_k with c_k = E_1(r_{k-1}), one extraction for all
# the channels, so K is common to them. Before extraction k the decomposition
# ends with "terminalResidue" when every projection of r_{k-1} on the rows of
# directions has at most one extremum, with "modeLimit" when k - 1 equals
# maxModes, and with "insufficientExtrema" when E_1 is undefined on r_{k-1}
# because a direction has no envelope. modes is samples x modes x channels and
# the residue samples x channels. BEMD stops before adding an incomplete
# candidate; MEMD keeps its existing policy through acceptIncomplete.
.decomposeVector <- function(signal, directions, envelope, sifting, acceptIncomplete) {
    Residue <- signal
    Modes <- list()
    repeat {
        Iterations <- 0L
        Direction <- NA_integer_
        Counts <- apply(Residue %*% t(directions), 2L, function(projection) {
            Extrema <- .findExtrema(projection)
            nrow(Extrema$maxima) + nrow(Extrema$minima)
        })
        if (all(Counts <= 1L)) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == sifting$maxModes) {
            StopReason <- "modeLimit"
            break
        }
        Extraction <- .extractVectorIMF(Residue, directions = directions, envelope = envelope, sifting = sifting)
        if (is.null(Extraction$mode)) {
            StopReason <- "insufficientExtrema"
            Direction <- Extraction$direction
            break
        }
        if (!acceptIncomplete && identical(Extraction$stopReason, "insufficientExtrema")) {
            Iterations <- Extraction$iterations
            Direction <- Extraction$direction
            StopReason <- Extraction$stopReason
            break
        }
        Modes[[length(Modes) + 1L]] <- Extraction
        Residue <- Residue - Extraction$mode
    }
    OUT <- list(
        modes = aperm(array(as.double(unlist(lapply(Modes, `[[`, "mode"), use.names = FALSE)),
                            dim = c(dim(signal), length(Modes))), c(1L, 3L, 2L)),
        residue = Residue,
        iterations = vapply(Modes, `[[`, integer(1L), "iterations"),
        converged = vapply(Modes, `[[`, logical(1L), "converged"),
        stopReason = StopReason
    )
    if (!acceptIncomplete) {
        OUT$termination <- data.frame(mode = length(Modes) + 1L, direction = Direction,
                                      iterations = Iterations, stopReason = StopReason)
    }
    OUT
}
