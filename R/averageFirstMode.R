# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Average first mode <E_1(r + p^(m))>_M of eq-completas-promedio
# (ceemdan.model:23-27, :167-187) over the columns p^(m) of perturbations:
# noise of an ensemble stage or mask columns. sifting is the configuration
# already completed and validated by the public entry. When E_1 is undefined
# for a perturbed input, that column contributes zero and stays in the
# divisor M, so its local mean is the perturbed input, as in ceemdan.m:88-104
# of Colominas (revision af9a0f83). iterations sums the subtractions of the M
# extractions; converged is TRUE when all M accepted, FALSE when one ended
# without accepting or could not begin, and NA when all completed fixed
# sifting.
.averageFirstMode <- function(residue, perturbations, sifting) {
    Sum <- numeric(length(residue))
    Iterations <- 0L
    Converged <- logical(ncol(perturbations))
    Termination <- list()
    for (m in seq_len(ncol(perturbations))) {
        Extraction <- .extractIMF(residue + perturbations[, m], sifting = sifting)
        if (is.null(Extraction)) {
            Extrema <- .findExtrema(residue + perturbations[, m])
            Termination[[length(Termination) + 1L]] <- data.frame(
                phase = m, iterations = 0L,
                stopReason = if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L)
                    "terminalResidue" else "insufficientExtrema")
            next
        }
        if (isFALSE(Extraction$converged)) {
            Termination[[length(Termination) + 1L]] <- data.frame(
                phase = m, iterations = Extraction$iterations, stopReason = Extraction$stopReason)
        }
        Sum <- Sum + Extraction$mode
        Iterations <- Iterations + Extraction$iterations
        Converged[m] <- Extraction$converged
    }
    list(mode = Sum / ncol(perturbations), iterations = Iterations, converged = all(Converged),
         termination = if (length(Termination)) do.call(rbind, Termination) else
             data.frame(phase = integer(), iterations = integer(), stopReason = character()))
}
