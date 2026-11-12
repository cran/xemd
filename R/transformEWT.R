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

# Mirror extension of the EWT record and its DFT
# (descomposicion.model.es.qmd:117-121), with the convention of ewtpy EWT1D:
# for h = ceiling(N / 2), the samples x_{h-1}, ..., x_1 precede the record and
# x_{N-1}, ..., x_{N-h} follow it, so T = N + 2h - 1. observed holds the
# positions h, ..., h + N - 1 of the record inside the extension. The
# convention is EWT's own, distinct from the VMD extension of
# .transformReflection.
.transformEWT <- function(signal) {
    n <- length(signal)
    h <- ceiling(n / 2)
    list(spectrum = stats::fft(signal[c(rev(seq_len(h - 1)), seq_len(n),
                                        (n - 1):(n - h))]),
         observed = h - 1 + seq_len(n))
}
