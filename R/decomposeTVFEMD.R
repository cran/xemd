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

# Time varying filter based empirical mode decomposition (tvf.model:125-353)
# by successive extraction: r_0 is the signal; before extraction k the
# decomposition ends with "terminalResidue" when r_{k-1} has fewer than four
# extrema, maxima plus minima (pytvfemd.py:69), and with "modeLimit" when k - 1
# equals maxModes, Inf by default; otherwise c_k is the component of r_{k-1}
# (.extractTVFMode) and r_k = r_{k-1} - c_k. The components satisfy the local
# narrow-band criterion of TVF-EMD, not the two conditions of an IMF
# (tvf.model:131). Positions are sample indices: bandwidthThreshold xi and
# intermittencyThreshold rho are dimensionless and the cutoff is in rad per
# sample. The decomposition runs on the bare samples, so
# no attribute of signal (names, tsp, class) reaches the outputs.
decomposeTVFEMD <- function(signal, bandwidthThreshold, splineOrder, intermittencyThreshold,
                            maxSift, maxModes = Inf) {
    isPositive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
    isCount <- function(x) isPositive(x) && x >= 1 && x == floor(x)
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal))) {
        stop("signal must be a numeric vector of finite real values without dim")
    }
    if (!isPositive(bandwidthThreshold)) {
        stop("bandwidthThreshold must be a finite number greater than zero")
    }
    if (!isCount(splineOrder)) stop("splineOrder must be a whole number of at least one")
    if (!isPositive(intermittencyThreshold)) {
        stop("intermittencyThreshold must be a finite number greater than zero")
    }
    if (!isCount(maxSift)) stop("maxSift must be a whole number of at least one")
    if (!(is.numeric(maxModes) && length(maxModes) == 1L && !is.na(maxModes) &&
          maxModes >= 1 && maxModes == floor(maxModes))) {
        stop("maxModes must be a whole number of at least one or Inf")
    }
    signal <- as.double(signal)
    Residue <- signal
    Modes <- list()
    repeat {
        Extrema <- .findExtrema(Residue)
        if (nrow(Extrema$maxima) + nrow(Extrema$minima) < 4L) {
            StopReason <- "terminalResidue"
            break
        }
        if (length(Modes) == maxModes) {
            StopReason <- "modeLimit"
            break
        }
        Extraction <- .extractTVFMode(Residue, bandwidthThreshold = bandwidthThreshold,
                                      splineOrder = splineOrder,
                                      intermittencyThreshold = intermittencyThreshold,
                                      maxSift = maxSift, first = length(Modes) == 0L)
        Modes[[length(Modes) + 1L]] <- Extraction
        Residue <- Residue - Extraction$mode
    }
    list(
        modes = matrix(as.double(unlist(lapply(Modes, `[[`, "mode"), use.names = FALSE)),
                       nrow = length(signal), ncol = length(Modes)),
        residue = Residue,
        iterations = vapply(Modes, `[[`, integer(1L), "iterations"),
        converged = vapply(Modes, `[[`, logical(1L), "converged"),
        stopReason = StopReason
    )
}
