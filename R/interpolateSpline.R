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
# Interpolant of envelope nodes (emd.model:35-42). With four or more nodes it
# is the cubic spline closed by the not-a-knot conditions of interp1 'spline'
# in emd.m, through pracma::interp1; with three nodes the parabola through
# them and with two the line. The nodes must be strictly increasing and cover
# every evaluation position: nothing is extrapolated. An evaluation position
# that is a node returns exactly the node value, as the nodal conditions state
# (emd.model:29-31), instead of the rounded evaluation; the result has no names.
.interpolateSpline <- function(position, value, at) {
    k <- length(position)
    if (k < 2L || any(diff(position) <= 0) || any(at < position[1L] | at > position[k])) {
        stop("spline nodes must be at least two, strictly increasing and cover every evaluation position")
    }
    if (k == 2L) {
        OUT <- value[1L] + (at - position[1L]) * (value[2L] - value[1L]) / (position[2L] - position[1L])
    }
    if (k == 3L) {
        Slopes <- diff(value) / diff(position)
        Curvature <- (Slopes[2L] - Slopes[1L]) / (position[3L] - position[1L])
        OUT <- value[1L] + (at - position[1L]) * (Slopes[1L] + (at - position[2L]) * Curvature)
    }
    if (k >= 4L) OUT <- pracma::interp1(x = position, y = value, xi = at, method = "spline")
    IDX <- match(at, position)
    OUT[!is.na(IDX)] <- value[IDX[!is.na(IDX)]]
    unname(OUT)
}
