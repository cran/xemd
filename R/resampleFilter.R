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

# Adaptation of the prototype filter of 2 m + 1 weights to the support
# 2 halfWidth + 1 by integrating its weights over the rescaled intervals
# (lmd.model:424), with the discretization of get_mask_v1_1 of fifpy. To
# reduce it, bin i of 2 l + 1 spans s = (i - 1) n / (2 l + 1) + 1 to
# e = i n / (2 l + 1) over the n weights y and takes the weights ceiling(s),
# ..., floor(e), none when ceiling(s) > floor(e), plus the fractions
# (ceiling(s) - s) y[ceiling(s)] and (e - floor(e)) y[floor(e)]; the bins are
# not renormalized, as in fifpy, so their mass may differ slightly from that
# of the prototype. To enlarge it, the right half y[m + 1], ..., y[n] is
# interpolated linearly at l + 1 equally spaced points, mirrored, and
# normalized to unit absolute sum (eq-if-normalizacion, lmd.model:226-230).
.resampleFilter <- function(filter, halfWidth) {
    n <- length(filter)
    m <- (n - 1L) %/% 2L
    Size <- 2L * halfWidth + 1L
    if (halfWidth <= m) {
        return(vapply(seq_len(Size), function(i) {
            s <- (i - 1) * n / Size + 1
            e <- i * n / Size
            Inner <- if (ceiling(s) <= floor(e)) sum(filter[ceiling(s):floor(e)]) else 0
            Inner + (ceiling(s) - s) * filter[ceiling(s)] + (e - floor(e)) * filter[floor(e)]
        }, numeric(1L)))
    }
    Half <- stats::approx(x = 0:m, y = filter[(m + 1L):n], xout = seq(0, m, length.out = halfWidth + 1L))$y
    Weights <- c(rev(Half[-1L]), Half)
    Weights / sum(abs(Weights))
}
