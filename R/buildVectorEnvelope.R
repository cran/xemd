# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Envelope e_d of a vector signal at the maxima of one of its projections
# (memd.model:24-30, :88-94): the nodes are the maxima of the projection,
# extended by .extendExtrema, and the value of a node is the vector sample at
# the coordinate it comes from (column source). A plateau of even length puts
# that coordinate between two samples: each channel is then interpolated
# linearly between them, the vector value at the exact centre of the plateau.
# Every channel is interpolated through the same nodes by .interpolateSpline.
# NULL when the boundary rule cannot extend the extrema of the projection. The
# rule must be reflection, whose nodes all come from samples; the border
# extrapolation compares scalar values (eemd.model:283-293) and gives no
# vector value.
.buildVectorEnvelope <- function(signal, projection, extrema, boundary) {
    Nodes <- .extendExtrema(projection, extrema = extrema, boundary = boundary)
    if (is.null(Nodes)) return(NULL)
    Source <- Nodes$upper[, "source"]
    Value <- (signal[floor(Source), , drop = FALSE] + signal[ceiling(Source), , drop = FALSE]) / 2
    apply(Value, 2L, function(value) {
        .interpolateSpline(position = Nodes$upper[, "position"], value = value, at = seq_len(nrow(signal)))
    })
}
