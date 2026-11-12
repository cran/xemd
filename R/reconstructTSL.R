# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Reconstruction of each series of decomposition, the list of tables of
# decomposeTSL, from the modes that selection retains: their sum and, when
# includeResidue, the remainder of the series, residue or fitResidual; a
# series without retained modes sums zero, and a method without remainder
# (EWT) adds none. With every mode and the remainder it is the reconstruction
# identity of the decomposition (eq-marco-reconstruccion), which in EEMD with
# the divisor "realizations" is the signal plus meanNoise
# (eq-eemd-reconstruccion-variable). selection is a table with the keys and
# mode of every mode of the decomposition, once, and retained, TRUE or FALSE,
# as the selection of selectModesTSL. includeResidue = TRUE, the complete
# reconstruction, is its only default. Returns the data.table reconstruction
# with the keys of the decomposition, t and s, on the grid and in the order of
# signal. The inputs are not modified.
reconstructTSL <- function(decomposition, selection, includeResidue = TRUE) {
    Keys <- .checkDecomposition(decomposition = decomposition, columns = character())
    if (!(isTRUE(includeResidue) || isFALSE(includeResidue))) stop("includeResidue must be TRUE or FALSE")
    Modes <- unique(decomposition$modes[, c(Keys, "mode"), with = FALSE])
    if (!(data.table::is.data.table(selection) && all(c(Keys, "mode", "retained") %in% names(selection)) &&
          nrow(selection) == nrow(Modes) && !anyDuplicated(selection, by = c(Keys, "mode")) &&
          !anyNA(selection[Modes, on = c(Keys, "mode"), which = TRUE]))) {
        stop("selection must name every mode of decomposition once, as the selection of selectModesTSL")
    }
    if (!(is.logical(selection$retained) && !anyNA(selection$retained))) {
        stop("selection$retained must be TRUE or FALSE")
    }
    Grid <- c(Keys, "t")
    OUT <- decomposition$signal[, Grid, with = FALSE]
    Retained <- selection$retained[selection[decomposition$modes, on = c(Keys, "mode"), which = TRUE]]
    IDX <- factor(OUT[decomposition$modes, on = Grid, which = TRUE][Retained], levels = seq_len(nrow(OUT)))
    s <- vapply(split(decomposition$modes$s[Retained], IDX), sum, numeric(1L))
    Remainder <- intersect(c("residue", "fitResidual"), names(decomposition))
    if (includeResidue && length(Remainder)) {
        IDX <- OUT[decomposition[[Remainder]], on = Grid, which = TRUE]
        s[IDX] <- s[IDX] + decomposition[[Remainder]]$s
    }
    set(OUT, j = "s", value = unname(s))
    list(reconstruction = OUT)
}
