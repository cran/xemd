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

# Half width l of the filter of one IF extraction, in samples, proportional to
# N / E for the E extrema of the residue and fixed during the extraction
# (lmd.model:232, 244), by the rule 'ave' of fifpy get_mask_length:
# l = 2 round(stretch N / E). From the second extraction, a value that does
# not exceed the previous half width becomes ceiling(1.1 previous), with the
# product formed in double precision as in fifpy. previous is NULL in the
# first extraction.
.selectHalfWidth <- function(extremaCount, sampleCount, stretch, previous) {
    HalfWidth <- 2 * round(sampleCount / extremaCount * stretch)
    if (!is.null(previous) && HalfWidth <= previous) HalfWidth <- ceiling(previous * 1.1)
    as.integer(HalfWidth)
}
