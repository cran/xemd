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

# Empirical filter bank of EWT (ewt.sst.model.es.qmd:21-91) on the folded grid
# |omega_nu| = 2 pi min(nu, T - nu) / T rad/sample of a DFT of length
# T = size, one column per band: H_0 of eq-ewt-filtro-aproximacion and H_j of
# eq-ewt-filtro-banda between consecutive boundaries, the last band closed at
# Omega_K = pi (eq-ewt-fronteras), with the polynomial of eq-ewt-beta, the
# coordinate of eq-ewt-coordenada-transicion and tau_j = gamma Omega_j. gamma
# is the bound of eq-ewt-gamma over Omega_0 = 0, ..., Omega_K = pi times
# 1 - 1 / T, ewtpy's strictly smaller value, so that neighbouring transitions
# do not overlap (eq-ewt-separacion-transiciones) and the squares of the
# filters add up to one up to (1 - gamma) pi (eq-ewt-particion-unidad).
# Beyond, the falling transition of the last band stops at pi, as in ewtpy
# (ewt.sst.model.es.qmd:67). Without boundaries the bank is a single band.
.buildFilterBank <- function(boundaries, size) {
    if (length(boundaries) == 0L) return(matrix(1, nrow = size, ncol = 1L))
    nu <- seq_len(size) - 1
    AngularFrequency <- 2 * pi * pmin(nu, size - nu) / size
    Edges <- c(0, boundaries, pi)
    Gamma <- (1 - 1 / size) * min(diff(Edges) / (Edges[-1L] + Edges[-length(Edges)]))
    transitionAngle <- function(frequency, boundary) {
        u <- (frequency - (1 - Gamma) * boundary) / (2 * Gamma * boundary)
        pi * u^4 * (35 - 84 * u + 70 * u^2 - 20 * u^3) / 2
    }
    Bank <- matrix(0, nrow = size, ncol = length(Edges) - 1L)
    for (k in seq_len(ncol(Bank))) {
        Lower <- Edges[k]
        Upper <- Edges[k + 1L]
        MASK <- AngularFrequency >= (1 + Gamma) * Lower &
            AngularFrequency <= (1 - Gamma) * Upper
        Bank[MASK, k] <- 1
        MASK <- AngularFrequency > (1 - Gamma) * Upper &
            AngularFrequency <= (1 + Gamma) * Upper
        Bank[MASK, k] <- cos(transitionAngle(AngularFrequency[MASK], boundary = Upper))
        if (Lower > 0) {
            MASK <- AngularFrequency >= (1 - Gamma) * Lower &
                AngularFrequency < (1 + Gamma) * Lower
            Bank[MASK, k] <- sin(transitionAngle(AngularFrequency[MASK], boundary = Lower))
        }
    }
    Bank
}
