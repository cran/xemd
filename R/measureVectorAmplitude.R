# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Amplitude of the envelopes of one state of the vector sifting at every
# sample, with the Euclidean norm over the channels (Rehman and Mandic 2011,
# section II). envelopes holds one samples x channels matrix per direction, in
# the order of the rows of the directions, and mean is their mean.
# "globalMean" is (1 / D) sum_d |e_d - m|, the distance of the D envelopes to
# their mean; "envelopePairs" is (1 / D) sum_{d <= D/2} |e_d - e_{d + D/2}|,
# the half-distances between the envelopes of antipodal directions averaged
# over the D / 2 pairs, as in memd.m of Rehman and Mandic (lines 237, 245);
# it needs the second half of the directions to be the antipodes of the first
# in the same order.
.measureVectorAmplitude <- function(envelopes, mean, amplitudeRule) {
    D <- length(envelopes)
    if (amplitudeRule == "globalMean") {
        return(Reduce(`+`, lapply(envelopes, function(e) sqrt(rowSums((e - mean)^2)))) / D)
    }
    if (amplitudeRule == "envelopePairs") {
        IDX <- seq_len(D / 2)
        return(Reduce(`+`, Map(function(e, antipode) sqrt(rowSums((e - antipode)^2)),
                               envelopes[IDX], envelopes[IDX + D / 2])) / D)
    }
    stop("unrecognized amplitudeRule: ", amplitudeRule)
}
