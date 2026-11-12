# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-FileCopyrightText: 2022 Stefano Bianchi
# SPDX-FileCopyrightText: 2022 Alessandro Longo
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from pytvfemd 0.3.14 (Stefano Bianchi and Alessandro Longo), GPL-3;
# modified for xemd on 2026-10-11.
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2022-  Stefano Bianchi, Alessandro Longo

# Least-squares B-spline approximation of order n (degree n - 1) with the
# given increasing knots k_0 < ... < k_m (eq-tvf-spline, eq-tvf-ajuste-spline,
# tvf.model:135-148), in the periodic spline space of period p = k_m - k_0 of
# pytvfemd 0.3.14 (pytvfemd.py:476-631): the knots continue on each side with
# n - 1 knots whose spacings repeat periodically those of the other end, the
# coefficients of the basis functions that differ by m are equal, and every
# sample position t enters as ((t - k_0) mod p) + k_0. The sum over the samples
# is then the sum of eq-tvf-ajuste-spline over one period of the periodic
# extension of signal, and the approximation is evaluated at the same folded
# positions. Solved by QR; returns the fitted values at the samples.
.approximateSpline <- function(signal, knots, order) {
    Pieces <- length(knots) - 1L
    Spacing <- diff(knots)
    Cycle <- (seq_len(order - 1L) - 1L) %% Pieces + 1L
    Augmented <- c(rev(knots[1L] - cumsum(rev(Spacing)[Cycle])), knots,
                   knots[Pieces + 1L] + cumsum(Spacing[Cycle]))
    Position <- (seq_along(signal) - knots[1L]) %% (knots[Pieces + 1L] - knots[1L]) + knots[1L]
    Design <- splines::splineDesign(Augmented, Position, ord = order)
    Basis <- t(rowsum(t(Design), (seq_len(ncol(Design)) - 1L) %% Pieces))
    as.vector(qr.fitted(qr(Basis), signal))
}
