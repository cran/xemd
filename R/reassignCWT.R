# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979.
# MATLAB original: Eugene Brevdo; R port: Dongik Jang, Hee-Seok Oh and Donghoh Kim.
# Synchrosqueezing toolbox attribution also names Hau-Tieng Wu.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-06 and 2026-10-11.

# Synchrosqueezing of the CWT, eq-sst-frecuencia-compleja and
# eq-sst-reasignacion (ewt.sst.model.es.qmd:140-159), discretized as SynchWave
# (synchrosqueezing.implementation.es.qmd:51). Only the coefficients with
# |W| > threshold, the strict domain of eq-sst-reasignacion, enter: their
# frequency Re Omega_f / (2 pi) = Im(dW / W) / (2 pi), in Hz because the
# derivative is in rad/s, goes to the nearest of J logarithmic cells from
# 1 / ((N - 1) samplingInterval) to 1 / (2 samplingInterval), J the number of
# scales, where the coefficient adds W a^(-1/2), a = 2^(j / voiceCount); the
# sums are divided by voiceCount. Estimates not positive or outside the cells
# are not accumulated, and a CWT with sqrt(sum |W|^2 / 2) below the machine
# epsilon gives a zero transform. Result: transform, samples by cells, and the
# cell frequencies in Hz.
.reassignCWT <- function(coefficients, derivative, samplingInterval, voiceCount,
                         threshold) {
    n <- nrow(coefficients)
    ScaleCount <- ncol(coefficients)
    Limits <- log2(1 / (c(n - 1, 2) * samplingInterval))
    Frequency <- 2^seq(Limits[1L], Limits[2L], length.out = ScaleCount)
    Transform <- matrix(0i, nrow = n, ncol = ScaleCount)
    if (sqrt(sum(Mod(coefficients)^2) / 2) < .Machine$double.eps) {
        return(list(transform = Transform, frequency = Frequency))
    }
    for (j in seq_len(ScaleCount)) {
        Rows <- which(Mod(coefficients[, j]) > threshold)
        Estimate <- Im(derivative[Rows, j] / coefficients[Rows, j]) / (2 * pi)
        MASK <- Estimate > 0
        Rows <- Rows[MASK]
        Cell <- floor(0.5 + (ScaleCount - 1) * (log2(Estimate[MASK]) - Limits[1L]) /
                          (Limits[2L] - Limits[1L])) + 1
        MASK <- Cell >= 1 & Cell <= ScaleCount
        IDX <- cbind(Rows[MASK], Cell[MASK])
        Transform[IDX] <- Transform[IDX] +
            coefficients[Rows[MASK], j] * 2^(-j / (2 * voiceCount))
    }
    list(transform = Transform / voiceCount, frequency = Frequency)
}
