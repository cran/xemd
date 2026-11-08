# Adapted from pyewt, revision 9d60a340f34e92356ed4f7a6a8b7044f524455e3.
# SPDX-License-Identifier: MIT
# MIT License
# 
# Copyright (c) 2024 Jerome Gilles
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

# Components of EWT on the record (ewt.sst.model.es.qmd:103-115): each column
# of coefficients is filtered by the dual filter H_j / S, S = sum_i H_i^2, as
# the inverse iewt1d of pyewt does, and .restoreObserved inverts the DFT and
# keeps the record. With the real bank, column j is IDFT(X H_j^2 / S), which
# is eq-ewt-componente where S = 1, and the columns add up to the record. A
# frequency where every filter vanishes leaves the dual bank undefined.
.synthesizeBands <- function(coefficients, filters, observed) {
    Power <- rowSums(filters^2)
    if (any(Power == 0)) {
        stop("the filter bank leaves frequencies uncovered: its dual is undefined")
    }
    .restoreObserved(spectrum = stats::mvfft(coefficients) * filters / Power,
                     observed = observed)
}
