# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from LMD 1.2.1, revision ec11b721ec7b05f0137eeb67a619b28f62ef85a2.
# Author: Shubhra Prakash. Source license: Apache License (>= 2).
# Modified for xemd on 2026-10-11.

# One product function of LMD by multiplicative demodulation (eq-lmd-centrado,
# eq-lmd-normalizacion, eq-lmd-envolvente-acumulada, eq-lmd-producto;
# lmd.model:53-77). From s_0 = signal, iteration q first counts the extrema of
# s_{q-1}: with at most one the extraction cannot continue and ends with
# Q = q - 1 factors and converged FALSE. Otherwise the local levels of s_{q-1}
# smoothed into m_q and a_q give s_q = (s_{q-1} - m_q) / a_q; a magnitude that
# is not positive leaves s_q undefined (lmd.model:65) and raises the condition
# xemdZeroEnvelope. When mean |1 - a_q| <= envelopeTolerance or
# mean |s_{q-1} - s_q| <= changeTolerance the extraction ends without that
# normalization, with Q = q - 1 and converged TRUE, in the order of
# extract_product_function of LMD 1.2.1; after maxIterations normalizations
# Q = maxIterations and converged FALSE. The product function is
# a_1 ... a_Q s_Q and iterations is Q.
.extractProduct <- function(signal, maxIterations, envelopeTolerance, changeTolerance,
                            maxSmoothPasses) {
    N <- length(signal)
    Current <- signal
    Envelopes <- list()
    Converged <- FALSE
    for (q in seq_len(maxIterations)) {
        Extrema <- .findExtrema(Current)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) <= 1L) break
        Levels <- .buildLocalLevels(Current, extrema = Extrema)
        Mean <- .smoothLocalLevels(Levels$mean, halfWidth = Levels$halfWidth,
                                   maxSmoothPasses = maxSmoothPasses)
        Envelope <- .smoothLocalLevels(Levels$amplitude, halfWidth = Levels$halfWidth,
                                       maxSmoothPasses = maxSmoothPasses)
        if (any(Envelope <= 0)) {
            stop(errorCondition("the local magnitude is not positive, so the normalization is undefined",
                                class = c("xemdZeroEnvelope", "moduleFailure")))
        }
        Candidate <- (Current - Mean) / Envelope
        if (sum(abs(1 - Envelope)) / N <= envelopeTolerance ||
            sum(abs(Current - Candidate)) / N <= changeTolerance) {
            Converged <- TRUE
            break
        }
        Envelopes[[q]] <- Envelope
        Current <- Candidate
    }
    list(mode = Reduce(`*`, Envelopes, Current), iterations = length(Envelopes), converged = Converged)
}
