# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Successive scalar extraction with validated inputs. EMD keeps its existing
# incomplete candidates; EEMD and CEEMD stop each realization before adding a
# candidate whose sifting lost its envelopes, leaving the residue at the input
# of that attempt. A candidate that reaches the iteration limit is kept with
# converged FALSE, as the reference libraries do. Only the strict consumers
# receive the attempted mode, work and cause in termination.
.decomposeScalar <- function(signal, sifting, acceptIncomplete) {
    Residue <- signal
    Modes <- list()
    repeat {
        Iterations <- 0L
        Extrema <- .findExtrema(Residue)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == sifting$maxModes) {
            StopReason <- "modeLimit"
            break
        }
        Extraction <- .extractIMF(Residue, sifting = sifting)
        if (is.null(Extraction)) {
            StopReason <- "insufficientExtrema"
            break
        }
        if (!acceptIncomplete && identical(Extraction$stopReason, "insufficientExtrema")) {
            Iterations <- Extraction$iterations
            StopReason <- Extraction$stopReason
            break
        }
        Modes[[length(Modes) + 1L]] <- Extraction
        Residue <- Residue - Extraction$mode
    }
    OUT <- list(
        modes = matrix(as.double(unlist(lapply(Modes, `[[`, "mode"), use.names = FALSE)),
                       nrow = length(signal), ncol = length(Modes)),
        residue = Residue,
        iterations = vapply(Modes, `[[`, integer(1L), "iterations"),
        converged = vapply(Modes, `[[`, logical(1L), "converged"),
        stopReason = StopReason
    )
    if (!acceptIncomplete) {
        OUT$termination <- data.frame(mode = length(Modes) + 1L, iterations = Iterations,
                                      stopReason = StopReason)
    }
    OUT
}
