# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-FileCopyrightText: 2022 Emanuele Papini
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from fifpy, revision 26e51e7d45e551950dea8be37ece447fc4340296.
# SPDX-License-Identifier: MIT
# MIT License
#
# Copyright (c) 2022 Emanuele Papini
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

# One IF extraction with the periodic extension of the record, evaluated as
# FIF (eq-if-potencia-espectral, eq-fif-diagonalizacion, eq-fif-extraccion;
# lmd.model:250-253, 428-442): the convolution is the circulant W whose first
# column holds w_j at lag j modulo N, the 2 l + 1 weights folded on the N
# samples when the filter is longer than the record, and its eigenvalues
# lambda are the discrete Fourier transform of that column. The spectrum of
# f^(q) is (1 - lambda)^q times the spectrum of signal; after each q the
# relative change of eq-if-parada-relativa on the full spectra is compared
# with tolerance, accepting with <= as fifpy; a null previous state leaves it
# undefined and the predicate false. Returns f^(Q) with Q the first accepted q
# (converged TRUE) or maxIterations (converged FALSE).
.iterateFluctuationPeriodic <- function(signal, filter, tolerance, maxIterations) {
    N <- length(signal)
    HalfWidth <- (length(filter) - 1L) %/% 2L
    Blocks <- ceiling(length(filter) / N)
    Folded <- rowSums(matrix(c(filter, numeric(Blocks * N - length(filter))), nrow = N))
    Multiplier <- 1 - stats::fft(Folded[(seq_len(N) - 1L + HalfWidth) %% N + 1L])
    Current <- stats::fft(signal)
    Iterations <- as.integer(maxIterations)
    Converged <- FALSE
    for (q in seq_len(maxIterations)) {
        Previous <- Current
        Current <- Multiplier * Previous
        if (isTRUE(.measureRelativeChange(previous = Previous, current = Current) <= tolerance)) {
            Iterations <- q
            Converged <- TRUE
            break
        }
    }
    list(mode = Re(stats::fft(Current, inverse = TRUE)) / N, iterations = Iterations,
         converged = Converged)
}
