# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Multivariate empirical mode decomposition (memd.model:74-156) of the C >= 2
# columns of signal: the vector sifting with the trajectory envelopes on the
# 2D directions of .buildDirections. sifting admits the rules "fixed",
# "consecutive" and "doubleThreshold", with reflection only, because the
# border extrapolation compares scalar values (eemd.model:283-293); under
# "doubleThreshold" the field amplitudeRule, "globalMean" or "envelopePairs",
# selects the amplitude of the ratio. Defaults of memd.m of Rehman and Mandic,
# the code of the authors of MEMD: 64 directions (lines 21, 563-565) and the
# criterion 'stop' (:27-32, :567-568), the double threshold, with the
# amplitude of that code, between the envelopes of each direction (:237,
# :245); .checkSifting completes the other fields. The decomposition runs on
# the bare samples, so no attribute of signal reaches the outputs.
decomposeMEMD <- function(signal, directionCount = 64, sifting) {
    if (!(is.matrix(signal) && is.numeric(signal) && ncol(signal) >= 2L && all(is.finite(signal)))) {
        stop("signal must be a numeric matrix of finite real values with two or more columns")
    }
    if (!(is.numeric(directionCount) && length(directionCount) == 1L && is.finite(directionCount) &&
          directionCount >= 1 && directionCount == round(directionCount))) {
        stop("directionCount must be a whole number of at least one")
    }
    if (is.list(sifting) && is.null(sifting[["stopRule"]])) sifting$stopRule <- "doubleThreshold"
    if (is.list(sifting) && identical(sifting[["stopRule"]], "doubleThreshold")) {
        AmplitudeRule <- if (is.null(sifting[["amplitudeRule"]])) "envelopePairs" else sifting[["amplitudeRule"]]
        sifting$amplitudeRule <- NULL
    }
    sifting <- .checkSifting(sifting)
    if (!(sifting$stopRule %in% c("fixed", "consecutive", "doubleThreshold"))) {
        stop("sifting$stopRule must be fixed, consecutive or doubleThreshold for MEMD")
    }
    if (sifting$boundary != "reflection") stop("sifting$boundary must be reflection for MEMD")
    if (sifting$stopRule == "doubleThreshold") {
        if (!(is.character(AmplitudeRule) && length(AmplitudeRule) == 1L &&
              AmplitudeRule %in% c("globalMean", "envelopePairs"))) {
            stop("sifting$amplitudeRule must be globalMean or envelopePairs")
        }
        sifting$amplitudeRule <- AmplitudeRule
    }
    .decomposeVector(matrix(as.vector(signal), nrow = nrow(signal)),
                     directions = .buildDirections(directionCount, channelCount = ncol(signal)),
                     envelope = "trajectory", sifting = sifting, acceptIncomplete = TRUE)
}
