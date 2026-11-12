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

# Local cut-off of TVF-EMD on a provisional signal (tvf.model:166-279), in
# rad per sample. A and omega are the amplitude and the frequency of the
# shared analytic signal, the frequency converted here from cycles per sample
# (samplingInterval 1) to rad per sample (eq-tvf-analitica).
# beta_1 and eta_1 interpolate A and the product omega A^2 at the minima of A,
# beta_2 and eta_2 at its maxima (eq-tvf-beta-minima, eq-tvf-beta-maxima,
# eq-tvf-eta-minima, eq-tvf-eta-maxima); the component amplitudes are
# a_1 = (beta_2 + beta_1) / 2 and a_2 = (beta_2 - beta_1) / 2
# (eq-tvf-amplitud-uno, eq-tvf-amplitud-dos) and the cutoff is
# (eta_2 - eta_1) / (4 a_1 a_2), NA where a_1 a_2 = 0 (eq-tvf-corte,
# tvf.model:279). The interpolation is the pchip of SciPy that pytvfemd 0.3.14
# applies (pytvfemd.py:299-310), through pracma::pchip, which extrapolates
# with its end pieces and needs three nodes; through two nodes it is the line.
# A node at a half-integer plateau centre takes the mean frequency of its two
# central samples. Returns NULL when A has fewer than two minima or two maxima
# (pytvfemd.py:288) or a node frequency is undefined (imf.model:124).
.estimateCutoff <- function(signal) {
    OUT <- computeInstantaneousAttributes(matrix(signal), samplingInterval = 1)
    Amplitude <- OUT$amplitude[, 1L]
    Frequency <- 2 * pi * OUT$frequency[, 1L]
    Extrema <- .findExtrema(Amplitude)
    if (nrow(Extrema$minima) < 2L || nrow(Extrema$maxima) < 2L) return(NULL)
    Product <- lapply(Extrema[c("minima", "maxima")], function(nodes) {
        Position <- nodes[, "position"]
        (Frequency[floor(Position)] + Frequency[ceiling(Position)]) / 2 * nodes[, "value"]^2
    })
    if (anyNA(unlist(Product, use.names = FALSE))) return(NULL)
    Sample <- seq_along(signal)
    interpolate <- function(nodes, value) {
        Position <- nodes[, "position"]
        if (length(Position) == 2L) {
            return(value[1L] + (Sample - Position[1L]) * diff(value) / diff(Position))
        }
        pracma::pchip(Position, value, Sample)
    }
    Beta1 <- interpolate(Extrema$minima, Extrema$minima[, "value"])
    Beta2 <- interpolate(Extrema$maxima, Extrema$maxima[, "value"])
    Eta1 <- interpolate(Extrema$minima, Product$minima)
    Eta2 <- interpolate(Extrema$maxima, Product$maxima)
    Components <- cbind(Beta2 + Beta1, Beta2 - Beta1, deparse.level = 0) / 2
    Denominator <- 4 * Components[, 1L] * Components[, 2L]
    Cutoff <- (Eta2 - Eta1) / Denominator
    Cutoff[Denominator == 0] <- NA_real_
    list(cutoff = Cutoff, amplitude = Components,
         eta = cbind(Eta1, Eta2, deparse.level = 0), analyticAmplitude = Amplitude)
}
