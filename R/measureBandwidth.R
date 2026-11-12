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

# Relative instantaneous bandwidth theta = B / mean frequency of TVF-EMD
# (eq-tvf-frecuencia-ponderada, eq-tvf-ancho-banda, eq-tvf-criterio-banda,
# tvf.model:305-337) from the estimate of .estimateCutoff, in rad per sample.
# The mean frequency is (eta_1 + eta_2) / (2 (a_1^2 + a_2^2)), the sum of the
# rows of eq-tvf-sistema-frecuencias (pytvfemd.py:326), and omega_1 - omega_2
# is the instantaneous frequency, through the shared analytic signal, of
# cos(phi_1 - phi_2) = (A^2 - a_1^2 - a_2^2) / (2 a_1 a_2) of
# eq-tvf-amplitud-compuesta restricted to [-1, 1] (pytvfemd.py:327-330): both
# avoid the separate frequencies, singular for equal amplitudes
# (tvf.model:279, :329). a_i' is the central difference per sample, each
# border repeating its neighbour (pytvfemd.py:332-335). theta is NA where the
# mean frequency is not positive or undefined, and at every sample when the
# cosine is undefined at one, because its analytic signal then is.
.measureBandwidth <- function(estimate) {
    Components <- estimate$amplitude
    n <- nrow(Components)
    Energy <- rowSums(Components^2)
    Product <- Components[, 1L] * Components[, 2L]
    Cosine <- (estimate$analyticAmplitude^2 - Energy) / (2 * Product)
    if (anyNA(Cosine)) return(rep(NA_real_, n))
    Cosine <- pmin(pmax(Cosine, -1), 1)
    Separation <- 2 * pi * computeInstantaneousAttributes(matrix(Cosine), samplingInterval = 1)$frequency[, 1L]
    Slope <- (Components[-c(1L, 2L), , drop = FALSE] - Components[-c(n - 1L, n), , drop = FALSE]) / 2
    Slope <- Slope[c(1L, seq_len(n - 2L), n - 2L), , drop = FALSE]
    Width <- sqrt(rowSums(Slope^2) / Energy + Product^2 * Separation^2 / Energy^2)
    Mean <- rowSums(estimate$eta) / (2 * Energy)
    OUT <- Width / Mean
    OUT[is.na(OUT) | is.na(Mean) | Mean <= 0] <- NA_real_
    OUT
}
