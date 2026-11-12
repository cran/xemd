# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-FileCopyrightText: 2019 Vinícius Carvalho
# SPDX-FileCopyrightText: 2019 Eduardo Mazoni
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from ewtpy, revision f9ff38bc5cd7859eea58c5b29099f8a503792b63.
# SPDX-License-Identifier: MIT
# Copyright (c) 2019 Vinícius Carvalho & Eduardo Mazoni
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

# Boundaries of the EWT partition (ewt.sst.model.es.qmd:9-17), in rad/sample,
# on the one-sided DFT magnitude s = |X_0|, ..., |X_{n-1}|, n = ceiling(N / 2).
# The maxima are those of .findExtrema, a plateau at its centre. The midpoint
# rule retains the bandCount - 1 largest maxima, or all of them when fewer
# exist, so that the number of bands adapts to K = M + 1
# (ewt.sst.model.es.qmd:17); with b_0 = 0 and b_1 < ... < b_{K-1} the 0-based
# bins of the retained maxima, Omega_j = pi (b_{j-1} + b_j) / N is the midpoint
# of the DFT frequencies 2 pi b / N (descomposicion.model.es.qmd:117-121,
# :135); ewtpy adds one bin. Equal values at the last retained maximum and the
# first discarded one leave the selection undefined.
#
# detection = NULL applies the midpoint rule to s. A detection list, completed
# and checked by decomposeEWT, applies the boundary detection of the EWT code
# of Gilles (supplement ewt, eq-deferred-ewt-filter and
# eq-deferred-ewt-threshold). With bandCount = NULL, the number of bands is
# the count of maxima of s with height above M_m + thresholdRatio (M_1 - M_m),
# M_1 and M_m the largest and smallest heights, and the midpoint rule follows.
# Otherwise the natural logarithm of s, which needs s > 0, precedes the
# regularization p = conv(s, h, 'same'), the central n values of the full
# convolution from position floor(L / 2) + 1, with the weights h_r = 1 / L of
# "average" or, for "gaussian", exp(-u_r^2 / (2 sigma^2)) normalized,
# u_r = r - (L - 1) / 2, r = 0, ..., L - 1, L = window; "none" keeps p = s.
# "locmax" applies the midpoint rule to p. "locmaxmin" and "locmaxminf" need
# bandCount - 1 positive maxima of p and take the largest, by frequency; on the
# interval from the first cell to the first of them and on those between
# consecutive ones, endpoints included, the boundary is the lowest minimum of
# p ("locmaxmin") or of s ("locmaxminf"), at the upper central position of a
# group of equal minima, and a minimum at position q gives 2 pi (q - 1) / N. An
# interval without minimum, or boundaries that do not increase, leave the
# partition undefined; no boundary is completed.
.segmentSpectrum <- function(signal, bandCount, detection = NULL) {
    n <- length(signal)
    Spectrum <- Mod(stats::fft(signal))[seq_len(ceiling(n / 2))]
    selectMaxima <- function(spectrum, count) {
        Maxima <- .findExtrema(spectrum)$maxima
        Maxima <- Maxima[order(Maxima[, "value"], decreasing = TRUE), , drop = FALSE]
        Count <- min(nrow(Maxima), count)
        if (Count > 0L && nrow(Maxima) > Count &&
            Maxima[Count + 1, "value"] == Maxima[Count, "value"]) {
            stop("tied spectral maxima leave the selection of ", count,
                 " maxima undefined")
        }
        sort(Maxima[seq_len(Count), "position"])
    }
    placeMidpoints <- function(positions) {
        Bins <- positions - 1
        pi * (c(0, utils::head(Bins, -1L)) + Bins) / n
    }
    if (is.null(detection)) return(placeMidpoints(selectMaxima(Spectrum, bandCount - 1)))
    if (is.null(bandCount)) {
        Heights <- .findExtrema(Spectrum)$maxima[, "value"]
        if (length(Heights) == 0L) stop("the spectrum has no maximum: zero bands")
        bandCount <- sum(Heights > min(Heights) +
                             detection$thresholdRatio * (max(Heights) - min(Heights)))
        if (bandCount == 0L) stop("no spectral maximum exceeds the relative threshold: zero bands")
        return(placeMidpoints(selectMaxima(Spectrum, bandCount - 1)))
    }
    if (detection$logarithm) {
        if (any(Spectrum <= 0)) stop("the logarithm needs a positive spectral magnitude")
        Spectrum <- log(Spectrum)
    }
    Regularized <- Spectrum
    if (detection$regularization != "none") {
        Window <- detection$window
        if (detection$regularization == "average") Weights <- rep(1 / Window, Window)
        if (detection$regularization == "gaussian") {
            Weights <- exp(-(seq_len(Window) - 1 - (Window - 1) / 2)^2 /
                               (2 * detection$sigma^2))
            Weights <- Weights / sum(Weights)
        }
        Padded <- c(numeric(Window), Spectrum, numeric(Window))
        Offset <- Window + floor(Window / 2) + 1 + seq_along(Spectrum)
        Regularized <- numeric(length(Spectrum))
        for (r in seq_len(Window)) Regularized <- Regularized + Weights[r] * Padded[Offset - r]
    }
    if (detection$detector == "locmax") return(placeMidpoints(selectMaxima(Regularized, bandCount - 1)))
    if (sum(.findExtrema(Regularized)$maxima[, "value"] > 0) < bandCount - 1) {
        stop("fewer than ", bandCount - 1, " positive maxima for the minima detector")
    }
    Edges <- c(1, selectMaxima(Regularized, bandCount - 1))
    if (detection$detector == "locmaxmin") Searched <- Regularized
    if (detection$detector == "locmaxminf") Searched <- Spectrum
    Minima <- .findExtrema(Searched)$minima
    Boundaries <- vapply(seq_len(length(Edges) - 1L), function(j) {
        Candidates <- Minima[Minima[, "position"] >= Edges[j] &
                                 Minima[, "position"] <= Edges[j + 1L], , drop = FALSE]
        if (nrow(Candidates) == 0L) {
            stop("no spectral minimum between positions ", Edges[j], " and ", Edges[j + 1L])
        }
        Lowest <- Candidates[Candidates[, "value"] == min(Candidates[, "value"]), "position"]
        Lowest[floor(length(Lowest) / 2) + 1]
    }, numeric(1L))
    if (any(diff(Boundaries) <= 0)) stop("the detected boundaries do not increase")
    2 * pi * (Boundaries - 1) / n
}
