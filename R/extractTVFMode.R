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

# One TVF-EMD component by sifting (tvf.model:339-347): h_0 is signal and, while
# the stopping rule fails, h_{s+1} = h_s - m_s (eq-tvf-sifting), m_s the
# periodic B-spline approximation of order splineOrder of h_s on the knots of
# its realigned cutoff. Every iteration extends h_s by symmetric reflection of
# round(N / 2) samples per side and restores the record afterwards
# (pytvfemd.py:75, :81-83, :162-163). The stopping rule is the aggregation of
# pytvfemd.py:111-126 in place of the local statement of tvf.model:337: with
# theta of .measureBandwidth on the samples round(0.05 N), ..., round(0.95 N)
# of the record, h_s is accepted when the mean of theta is below
# bandwidthThreshold (1 + (s + 1) / 4), from s = 1 in the first component and
# from s = 0 in the later ones, or when the samples with theta above
# bandwidthThreshold are fewer than 0.2 N; both run over the samples where
# theta is defined. converged is TRUE when the rule accepts and FALSE when
# maxSift subtractions are done without it or when the cutoff or the knots of
# h_s cannot be built; h_s is then the component. iterations counts them.
.extractTVFMode <- function(signal, bandwidthThreshold, splineOrder, intermittencyThreshold,
                            maxSift, first) {
    n <- length(signal)
    Padding <- round(n / 2)
    Record <- Padding + seq_len(n)
    Window <- round(0.05 * n):round(0.95 * n)
    Mode <- signal
    for (s in 0:maxSift) {
        Extended <- c(rev(Mode[seq_len(Padding) + 1L]), Mode, Mode[n - seq_len(Padding)])
        Estimate <- .estimateCutoff(Extended)
        if (is.null(Estimate)) return(list(mode = Mode, iterations = s, converged = FALSE))
        Bandwidth <- .measureBandwidth(Estimate)[Record][Window]
        if (((s >= 1L || !first) &&
             isTRUE(mean(Bandwidth, na.rm = TRUE) < bandwidthThreshold * (1 + (s + 1) / 4))) ||
            sum(Bandwidth > bandwidthThreshold, na.rm = TRUE) / n < 0.2) {
            return(list(mode = Mode, iterations = s, converged = TRUE))
        }
        if (s == maxSift) return(list(mode = Mode, iterations = s, converged = FALSE))
        Cutoff <- .realignCutoff(Estimate$cutoff, Extended, threshold = intermittencyThreshold)
        if (is.null(Cutoff)) return(list(mode = Mode, iterations = s, converged = FALSE))
        Knots <- .placeKnots(Cutoff)
        if (length(Knots) < 3L) return(list(mode = Mode, iterations = s, converged = FALSE))
        Mode <- Mode - .approximateSpline(Extended, knots = Knots, order = splineOrder)[Record]
    }
}
