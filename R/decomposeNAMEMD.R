# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Noise-assisted MEMD (memd.model:166-186): one MEMD of the signal augmented
# with the Q supplied noise channels, in the units of the signal, keeping the
# C coordinates of the signal, c_k = P_x c~_k. Projections, stops, iterations
# and termination are those of the augmented decomposition; the kernel draws
# no random number. directionCount and sifting are those of decomposeMEMD.
# The augmented matrix joins the bare samples, so that neither a time series
# method of cbind nor any attribute of the inputs reaches the outputs.
decomposeNAMEMD <- function(signal, noiseChannels, directionCount = 64, sifting) {
    if (!(is.matrix(signal) && is.numeric(signal) && ncol(signal) >= 1L && all(is.finite(signal)))) {
        stop("signal must be a numeric matrix of finite real values with one or more columns")
    }
    if (!(is.matrix(noiseChannels) && is.numeric(noiseChannels) && nrow(noiseChannels) == nrow(signal) &&
          ncol(noiseChannels) >= 1L && all(is.finite(noiseChannels)))) {
        stop("noiseChannels must be a numeric matrix of finite real values with the rows of signal ",
             "and one or more columns")
    }
    OUT <- decomposeMEMD(cbind(matrix(as.vector(signal), nrow = nrow(signal)),
                               matrix(as.vector(noiseChannels), nrow = nrow(noiseChannels))),
                         directionCount = directionCount, sifting = sifting)
    OUT$modes <- OUT$modes[, , seq_len(ncol(signal)), drop = FALSE]
    OUT$residue <- OUT$residue[, seq_len(ncol(signal)), drop = FALSE]
    OUT
}
