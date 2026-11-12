# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Selection of the modes of decomposition, the list of tables of decomposeTSL,
# with a declared criterion on metrics, the table metrics of measureModesTSL
# for the same decomposition. criterion is a list with rule, exactly the
# fields of that rule, and action:
# - rule "modes", field modes: a table with mode, whole numbers from 1, and
#   key columns of the decomposition; each row names the mode of that index of
#   every series whose keys it matches, so the indices go by group (by columns
#   and ID) or by series (with OCID), and each row names a mode of the
#   decomposition;
# - rule "energyFraction", field minimum: a mode meets it when its
#   energyFraction is not below minimum, a finite number not below zero;
# - rule "band", fields fMin and fMax in Hz, 0 <= fMin <= fMax: a mode meets
#   the closed band when fMin <= frequency <= fMax; a mode without frequency
#   (zero energy) is in no band. keep.Fo and remove.Fo of gmsp require
#   [Fo - BW, Fo + BW] inside the band (gmsp TS2IMF.R:218-239); this rule
#   compares frequency alone.
# action "keep" retains the modes that meet the criterion and "remove" the
# others. Returns the data.tables modes, the rows of the retained modes in the
# order of the decomposition, and selection, every mode with its keys, the
# criterion as declared text (for example "keep frequency in [2, 15] Hz") and
# retained. The inputs are not modified.
selectModesTSL <- function(decomposition, metrics, criterion) {
    Keys <- .checkDecomposition(decomposition = decomposition, columns = c("criterion", "retained"))
    Modes <- unique(decomposition$modes[, c(Keys, "mode"), with = FALSE])
    if (!(data.table::is.data.table(metrics) &&
          all(c(Keys, "mode", "frequency", "energyFraction") %in% names(metrics)) &&
          nrow(metrics) == nrow(Modes) && !anyDuplicated(metrics, by = c(Keys, "mode")) &&
          is.numeric(metrics$frequency) && is.numeric(metrics$energyFraction) &&
          !anyNA(metrics$energyFraction))) {
        stop("metrics must be the table metrics of measureModesTSL for the modes of decomposition")
    }
    IDX <- metrics[Modes, on = c(Keys, "mode"), which = TRUE]
    if (anyNA(IDX)) stop("metrics must be the table metrics of measureModesTSL for the modes of decomposition")
    Fields <- list(modes = "modes", energyFraction = "minimum", band = c("fMin", "fMax"))
    if (!(is.list(criterion) && is.character(criterion$rule) && length(criterion$rule) == 1L &&
          criterion$rule %in% names(Fields))) {
        stop("criterion$rule must be modes, energyFraction or band")
    }
    Rule <- criterion$rule
    if (!(length(names(criterion)) == length(criterion) && length(criterion) == length(Fields[[Rule]]) + 2L &&
          setequal(names(criterion), c("rule", Fields[[Rule]], "action")))) {
        stop("criterion of rule ", Rule, " must hold exactly the fields rule, ",
             paste(Fields[[Rule]], collapse = ", "), " and action")
    }
    Action <- criterion$action
    if (!(is.character(Action) && length(Action) == 1L && Action %in% c("keep", "remove"))) {
        stop("criterion$action must be keep or remove")
    }
    if (Rule == "modes") {
        Listed <- criterion$modes
        On <- c(setdiff(names(Listed), "mode"), "mode")
        if (!(is.data.frame(Listed) && "mode" %in% names(Listed) && !anyDuplicated(names(Listed)) &&
              all(On %in% c(Keys, "mode")))) {
            stop("criterion$modes must be a table of mode and key columns of decomposition")
        }
        Index <- Listed$mode
        if (!(is.numeric(Index) && all(is.finite(Index) & Index >= 1 & Index == round(Index)))) {
            stop("criterion$modes$mode must hold whole numbers from 1")
        }
        Listed <- data.table::as.data.table(Listed)
        set(Listed, j = "mode", value = as.integer(Index))
        if (anyNA(Modes[Listed, on = On, which = TRUE, mult = "first"])) {
            stop("criterion$modes names modes absent from decomposition")
        }
        Met <- !is.na(Listed[Modes, on = On, which = TRUE, mult = "first"])
        Label <- paste(Action, "listed modes")
    }
    if (Rule == "energyFraction") {
        Minimum <- criterion$minimum
        if (!(is.numeric(Minimum) && length(Minimum) == 1L && is.finite(Minimum) && Minimum >= 0)) {
            stop("criterion$minimum must be a finite number not below zero")
        }
        Met <- metrics$energyFraction[IDX] >= Minimum
        Label <- paste0(Action, " energyFraction >= ", Minimum)
    }
    if (Rule == "band") {
        Band <- c(criterion$fMin, criterion$fMax)
        if (!(is.numeric(Band) && length(Band) == 2L && all(is.finite(Band)) && Band[1L] >= 0 &&
              Band[1L] <= Band[2L])) {
            stop("criterion$fMin and criterion$fMax must be finite frequencies with 0 <= fMin <= fMax")
        }
        Frequency <- metrics$frequency[IDX]
        Met <- !is.na(Frequency) & Frequency >= Band[1L] & Frequency <= Band[2L]
        Label <- paste0(Action, " frequency in [", Band[1L], ", ", Band[2L], "] Hz")
    }
    Retained <- switch(Action, keep = Met, remove = !Met)
    Rows <- Modes[decomposition$modes, on = c(Keys, "mode"), which = TRUE]
    list(modes = decomposition$modes[Retained[Rows]],
         selection = cbind(Modes, data.table(criterion = Label, retained = Retained)))
}
