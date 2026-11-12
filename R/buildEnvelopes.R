# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Upper and lower envelopes e+ and e- at the positions 1, ..., N
# (emd.model:27-33): the nodes of .extendExtrema interpolated by
# .interpolateSpline. NULL when the boundary rule cannot extend the extrema.
# The envelopes pass through their nodes; they need not bound the signal.
.buildEnvelopes <- function(signal, extrema, boundary) {
    Nodes <- .extendExtrema(signal, extrema = extrema, boundary = boundary)
    if (is.null(Nodes)) return(NULL)
    lapply(Nodes, function(x) {
        .interpolateSpline(position = x[, "position"], value = x[, "value"], at = seq_along(signal))
    })
}
