# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Synchrosqueezed wavelet transform (ewt.sst.model.es.qmd:125-159) of a signal
# sampled every samplingInterval seconds: the Morlet CWT with voiceCount
# voices per octave and centre morletCenter, reassigned to frequency cells
# from the coefficients whose modulus exceeds threshold, in the units of the
# signal. Returns transform, complex samples by frequencies, and frequency in
# Hz; the representation recovers the signal only approximately.
transformSST <- function(signal, samplingInterval, voiceCount = 16,
                         morletCenter = 2 * pi, threshold) {
    stopifnot(
        "signal must be a numeric vector of at least four finite values" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 4L &&
            all(is.finite(signal)),
        "samplingInterval must be a finite number greater than zero" =
            is.numeric(samplingInterval) && length(samplingInterval) == 1L &&
            is.finite(samplingInterval) && samplingInterval > 0,
        "voiceCount must be a whole number not below one" =
            is.numeric(voiceCount) && length(voiceCount) == 1L &&
            is.finite(voiceCount) && voiceCount >= 1 &&
            voiceCount == floor(voiceCount),
        "morletCenter must be a finite number greater than zero" =
            is.numeric(morletCenter) && length(morletCenter) == 1L &&
            is.finite(morletCenter) && morletCenter > 0,
        "threshold must be a finite number greater than zero" =
            is.numeric(threshold) && length(threshold) == 1L &&
            is.finite(threshold) && threshold > 0)
    # The outputs inherit no attribute of signal (names, tsp, class).
    signal <- as.vector(signal)
    CWT <- .transformCWT(signal = signal, samplingInterval = samplingInterval,
                         voiceCount = voiceCount, morletCenter = morletCenter)
    .reassignCWT(coefficients = CWT$coefficients, derivative = CWT$derivative,
                 samplingInterval = samplingInterval, voiceCount = voiceCount,
                 threshold = threshold)
}
