# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Extrema and zero crossings of a sampled function (emd.model:17-25;
# imf.model:5-9). A maximal run of equal samples that touches neither sample 1
# nor sample N is a maximum (minimum) when both exterior neighbours are lower
# (higher); its position is the exact centre (i + j) / 2 and its value h_i, so
# an isolated extremum is the run of length one. A crossing is a sign change
# between neighbours or an interior run of zeros whose neighbours differ in
# sign; the counts N_e and N_z are the rows of maxima and minima and crossings.
.findExtrema <- function(signal) {
    Runs <- rle(unname(signal))
    Ends <- cumsum(Runs$lengths)
    IDX <- seq_along(Ends)[-c(1L, length(Ends))]
    Value <- Runs$values[IDX]
    Before <- Runs$values[IDX - 1L]
    After <- Runs$values[IDX + 1L]
    Position <- (Ends[IDX - 1L] + 1L + Ends[IDX]) / 2
    MASK.max <- Before < Value & After < Value
    MASK.min <- Before > Value & After > Value
    list(
        maxima = cbind(position = Position[MASK.max], value = Value[MASK.max]),
        minima = cbind(position = Position[MASK.min], value = Value[MASK.min]),
        crossings = sum(sign(signal[-1L]) * sign(signal[-length(signal)]) < 0) +
            sum(Value == 0 & sign(Before) * sign(After) < 0)
    )
}
