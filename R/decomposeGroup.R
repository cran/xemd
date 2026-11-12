# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Decomposition of one group of decomposeTSL, a record and an ID: series is
# the N x C matrix of its OCID in ascending order and samplingInterval the
# step of its grid in seconds. The operation derives the arguments that
# control does not carry: samplingInterval for the methods that declare it;
# for the ensembles, noise from realizationCount, one matrix of
# .generateNoise for each OCID in ascending order; for NA-MEMD, noiseChannels
# from noiseChannelCount, drawn by .generateNoise and scaled by noiseRatio
# times the reference deviation sqrt(mean(s_c^2)) of the sample deviations
# s_c of the OCID, without standardizing them. A scalar method decomposes each
# column; a coupled method decomposes the matrix, and its modes and its
# residue or fitResidual split by channel while its other fields, common to
# the channels, repeat. Returns the output of the method for each OCID.
.decomposeGroup <- function(series, samplingInterval, method, control) {
    Formals <- names(formals(method))
    Arguments <- control[intersect(names(control), Formals)]
    if ("samplingInterval" %in% Formals) Arguments$samplingInterval <- samplingInterval
    if ("noiseChannels" %in% Formals) {
        Arguments$noiseChannels <- control[["noiseRatio"]] * sqrt(mean(apply(series, 2L, stats::sd)^2)) *
            .generateNoise(nrow(series), control[["noiseChannelCount"]])
    }
    if (method %in% c("decomposeBEMD", "decomposeMEMD", "decomposeNAMEMD", "decomposeMVMD")) {
        AUX <- do.call(method, c(list(signal = series), Arguments))
        lapply(seq_len(ncol(series)), function(c) {
            Channel <- AUX
            Channel$modes <- matrix(AUX$modes[, , c], nrow = nrow(series))
            for (Field in intersect(c("residue", "fitResidual"), names(AUX))) {
                Channel[[Field]] <- AUX[[Field]][, c]
            }
            Channel
        })
    } else {
        lapply(seq_len(ncol(series)), function(c) {
            if ("noise" %in% Formals) {
                Arguments$noise <- .generateNoise(nrow(series), control[["realizationCount"]])
            }
            do.call(method, c(list(signal = series[, c]), Arguments))
        })
    }
}
