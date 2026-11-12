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

# Knots of the TVF-EMD approximation (eq-tvf-oscilacion-nudos,
# tvf.model:297-303): the positions of the maxima and the minima of
# q = cos(phi), phi the integral of the realigned cutoff in rad per sample
# from t_0 at the first sample, phi_1 = 0 and phi_{n+1} = phi_n + cutoff_n, as
# pytvfemd 0.3.14 does (pytvfemd.py:129-135). Increasing positions.
.placeKnots <- function(cutoff) {
    Extrema <- .findExtrema(cos(cumsum(c(0, cutoff[-length(cutoff)]))))
    sort(c(Extrema$maxima[, "position"], Extrema$minima[, "position"]))
}
