# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-FileCopyrightText: 2020 Marcelo A. Colominas
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from emd.m (G. Rilling, March 2007), distributed with CEEMDAN af9a0f83.
# SPDX-License-Identifier: MIT
# MIT License
# Copyright (c) 2020 Marcelo A. Colominas
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# Envelope nodes of each family of extrema (emd.model:228-240;
# eemd.model:283-295): the interior extrema of the family and the nodes that
# extend them to cover samples 1 to N, in strictly increasing position, as a
# matrix with exactly the columns position, value and source and no row names.
# source is the sample the node comes from: its own position for an interior
# extremum or a border node, tau for an extremum reflected to 2b - tau.
# Reflection keeps the selection of boundary_conditions in emd.m: two extrema
# per family and side, an axis on the first extremum or on the border sample,
# moved to the border when the reflected nodes do not cover it; the right
# border is the same rule in reversed time. Without three extrema it returns
# NULL. Extrapolation anchors each family at both borders with max{L+, h} and
# min{L-, h}, L+- the line through its two nearest extrema, or h with fewer
# than two: a border value need not be h, and source only names the border.
.extendExtrema <- function(signal, extrema, boundary) {
    n <- length(signal)
    buildNodes <- function(position, value, source) {
        matrix(c(position, value, source), ncol = 3L,
               dimnames = list(NULL, c("position", "value", "source")))
    }
    Maxima <- buildNodes(position = extrema$maxima[, "position"], value = extrema$maxima[, "value"],
                         source = extrema$maxima[, "position"])
    Minima <- buildNodes(position = extrema$minima[, "position"], value = extrema$minima[, "value"],
                         source = extrema$minima[, "position"])
    if (boundary == "reflection") {
        if (nrow(Maxima) + nrow(Minima) < 3L) return(NULL)
        reflect <- function(x, axis) {
            x <- x[rev(seq_len(nrow(x))), , drop = FALSE]
            x[, "position"] <- 2 * axis - x[, "position"]
            x
        }
        # Time reversal of the sample positions, its own inverse: it maps the
        # extrema to the left border and the reflected nodes back.
        reverse <- function(x) {
            x <- x[rev(seq_len(nrow(x))), , drop = FALSE]
            x[, c("position", "source")] <- n + 1 - x[, c("position", "source")]
            x
        }
        reflectLeft <- function(sample, maxima, minima) {
            Edge <- buildNodes(position = 1, value = sample, source = 1)
            pick <- function(x, rows) x[rows[rows <= nrow(x)], , drop = FALSE]
            MaximumFirst <- maxima[1L, "position"] < minima[1L, "position"]
            if (MaximumFirst && sample > minima[1L, "value"]) {
                Axis <- maxima[1L, "position"]
                Upper <- pick(maxima, 2:3)
                Lower <- pick(minima, 1:2)
            }
            if (MaximumFirst && sample <= minima[1L, "value"]) {
                Axis <- 1
                Upper <- pick(maxima, 1:2)
                Lower <- rbind(Edge, pick(minima, 1L))
            }
            if (!MaximumFirst && sample < maxima[1L, "value"]) {
                Axis <- minima[1L, "position"]
                Upper <- pick(maxima, 1:2)
                Lower <- pick(minima, 2:3)
            }
            if (!MaximumFirst && sample >= maxima[1L, "value"]) {
                Axis <- 1
                Upper <- rbind(Edge, pick(maxima, 1L))
                Lower <- pick(minima, 1:2)
            }
            if (nrow(Upper) == 0L || nrow(Lower) == 0L) {
                stop("reflection selected no extrema of a family with three or more extrema")
            }
            if (2 * Axis - max(Upper[, "position"]) > 1 || 2 * Axis - max(Lower[, "position"]) > 1) {
                if (Axis == 1) stop("reflection about the border sample does not cover the border")
                if (Axis == maxima[1L, "position"]) Upper <- pick(maxima, 1:2)
                if (Axis == minima[1L, "position"]) Lower <- pick(minima, 1:2)
                Axis <- 1
            }
            list(upper = reflect(Upper, axis = Axis), lower = reflect(Lower, axis = Axis))
        }
        Left <- reflectLeft(signal[[1L]], maxima = Maxima, minima = Minima)
        Right <- reflectLeft(signal[[n]], maxima = reverse(Maxima), minima = reverse(Minima))
        return(list(upper = rbind(Left$upper, Maxima, reverse(Right$upper)),
                    lower = rbind(Left$lower, Minima, reverse(Right$lower))))
    }
    if (boundary == "extrapolation") {
        anchor <- function(x, select) {
            k <- nrow(x)
            Ends <- c(signal[[1L]], signal[[n]])
            if (k >= 2L) {
                Ends <- select(Ends, c(
                    x[1L, "value"] + (1 - x[1L, "position"]) *
                        (x[2L, "value"] - x[1L, "value"]) / (x[2L, "position"] - x[1L, "position"]),
                    x[k, "value"] + (n - x[k, "position"]) *
                        (x[k, "value"] - x[k - 1L, "value"]) / (x[k, "position"] - x[k - 1L, "position"])
                ))
            }
            rbind(buildNodes(position = 1, value = Ends[1L], source = 1), x,
                  buildNodes(position = n, value = Ends[2L], source = n))
        }
        return(list(upper = anchor(Maxima, select = pmax), lower = anchor(Minima, select = pmin)))
    }
    stop("unrecognized boundary: ", boundary)
}
