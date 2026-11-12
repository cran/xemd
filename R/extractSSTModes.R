# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979:
# curve_ext_multi of the Synchrosqueezing Toolbox, which calls curve_ext
# (original author Jianfeng Lu, maintainer Eugene Brevdo); R port: Dongik Jang.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-15.

# Components of a synchrosqueezed transform along its ridges (supplement sst,
# eq-deferred-sst-band; curve_ext_multi and synsq_cwt_iw of the
# Synchrosqueezing Toolbox). transform is samples by frequency cells, as
# transformSST returns it, frequency its cell frequencies in Hz and
# morletCenter the Morlet centre that produced it. modeCount trajectories are
# extracted in succession by .traceSSTRidge with penalty; after each one the
# cells within exclusionWidth of its trajectory, clipped to the cells, take
# sqrt(.Machine$double.eps), so the next extraction renormalizes the energy.
# Component k is .invertSST of the original coefficients in the band of half
# width reconstructionWidth around trajectory k, clipped to the cells:
# Re sum_{j in B} T / C. The bands may overlap and no residue is added.
# Defaults: penalty = 1000 and exclusionWidth = 4, the values that
# curve_ext_multi assigns.
extractSSTModes <- function(transform, frequency, modeCount, penalty = 1000,
                            exclusionWidth = 4, reconstructionWidth, morletCenter) {
    isNumber <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    isCount <- function(x, lowest) isNumber(x) && x >= lowest && x == floor(x)
    stopifnot(
        "transform must be a numeric or complex matrix of finite values" =
            (is.numeric(transform) || is.complex(transform)) && is.matrix(transform) &&
            all(is.finite(transform)),
        "frequency must hold an increasing positive frequency per column of transform" =
            is.numeric(frequency) && is.null(dim(frequency)) &&
            length(frequency) == ncol(transform) && all(is.finite(frequency)) &&
            all(frequency > 0) && all(diff(frequency) > 0),
        "modeCount must be a whole number not below one" = isCount(modeCount, lowest = 1),
        "penalty must be a finite number not below zero" = isNumber(penalty) && penalty >= 0,
        "exclusionWidth must be a whole number not below zero" =
            isCount(exclusionWidth, lowest = 0),
        "reconstructionWidth must be a whole number not below zero" =
            isCount(reconstructionWidth, lowest = 0),
        "morletCenter must be a finite number greater than zero" =
            isNumber(morletCenter) && morletCenter > 0)
    inBand <- function(path, width) abs(outer(path, seq_len(ncol(transform)), "-")) <= width
    Work <- transform
    Modes <- matrix(0, nrow = nrow(transform), ncol = modeCount)
    for (k in seq_len(modeCount)) {
        Path <- .traceSSTRidge(transform = Work, frequency = frequency, penalty = penalty)
        Work[inBand(Path, width = exclusionWidth)] <- sqrt(.Machine$double.eps)
        Modes[, k] <- .invertSST(transform = transform * inBand(Path, width = reconstructionWidth),
                                 morletCenter = morletCenter)
    }
    Modes
}
