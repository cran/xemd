# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Masked decomposition by successive extraction (mascaras.model:106-180).
# r_0 is the signal; before stage k the decomposition ends, in this order, with
# "terminalResidue" when r_{k-1} has at most one extremum (emd.model:335), with
# "modeLimit" when k - 1 equals maxModes or the schedule is exhausted, and with
# "insufficientExtrema" when E_1 is undefined on r_{k-1} (emd.model:337).
# Otherwise schedule(k, r_{k-1}) gives the frequency nu_k in cycles per sample
# and the amplitude a_k of stage k; the masks are a_k cos(2 pi nu_k n + phase)
# with n = 0, ..., N - 1 (eq-mascaras-entradas), the first sample as phase
# origin as in get_next_imf_mask of emd 0.8.1 by Andrew Quinn
# (mascaras.implementation:211-213), the mode c_k is the average of
# E_1 over r_{k-1} plus each mask (eq-mascaras-upemd-primero,
# eq-mascaras-multinivel-modo) and r_k = r_{k-1} - c_k
# (eq-mascaras-multinivel-residuo). MaskedEMD rejects the whole stage if any
# phase is incomplete; itEMD retains its existing partial-stage policy.
.decomposeMasked <- function(signal, schedule, phases, sifting, acceptIncomplete) {
    Sample <- seq_along(signal) - 1
    Residue <- signal
    Modes <- list()
    Termination <- NULL
    repeat {
        Extrema <- .findExtrema(Residue)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) {
            StopReason <- "terminalResidue"
            break
        }
        Stage <- NULL
        if (length(Modes) < sifting$maxModes) Stage <- schedule(length(Modes) + 1L, Residue)
        if (is.null(Stage)) {
            StopReason <- "modeLimit"
            break
        }
        if (is.null(.buildEnvelopes(Residue, extrema = Extrema, boundary = sifting$boundary))) {
            StopReason <- "insufficientExtrema"
            break
        }
        Masks <- Stage$amplitude * cos(outer(2 * pi * Stage$cyclesPerSample * Sample, phases, "+"))
        Extraction <- .averageFirstMode(Residue, perturbations = Masks, sifting = sifting)
        if (!acceptIncomplete && any(Extraction$termination$stopReason != "iterationLimit")) {
            StopReason <- "incompleteExtraction"
            Termination <- data.frame(mode = length(Modes) + 1L, Extraction$termination)
            break
        }
        Modes[[length(Modes) + 1L]] <- c(Extraction, cyclesPerSample = Stage$cyclesPerSample)
        Residue <- Residue - Extraction$mode
    }
    OUT <- list(
        modes = matrix(as.double(unlist(lapply(Modes, `[[`, "mode"), use.names = FALSE)),
                       nrow = length(signal), ncol = length(Modes)),
        residue = Residue,
        iterations = vapply(Modes, `[[`, integer(1L), "iterations"),
        converged = vapply(Modes, `[[`, logical(1L), "converged"),
        stopReason = StopReason,
        cyclesPerSample = vapply(Modes, `[[`, double(1L), "cyclesPerSample")
    )
    if (!acceptIncomplete) {
        OUT$termination <- if (!is.null(Termination)) Termination else
            data.frame(mode = length(Modes) + 1L, phase = NA_integer_, iterations = 0L,
                       stopReason = StopReason)
    }
    OUT
}
