# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979.
# MATLAB original: Eugene Brevdo; R port: Dongik Jang, Hee-Seok Oh and Donghoh Kim.
# Synchrosqueezing toolbox attribution also names Hau-Tieng Wu.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-06 and 2026-10-11.

# Continuous wavelet transform eq-sst-cwt (ewt.sst.model.es.qmd:127-132) and
# its derivative in b, discretized as SynchWave with the Morlet wavelet
# (synchrosqueezing.implementation.es.qmd:3): the record is extended by
# symmetric padding, which repeats each edge sample, to
# T = 2^(1 + round(log2 N)) samples; the scales are a_j = 2^(j / voiceCount),
# j = 1, ..., (log2 T - 1) voiceCount; on the angular grid
# omega = 2 pi / T (0, ..., T / 2, -T / 2 + 1, ..., -1) rad/sample the filter
# of scale a is psi^(a omega) sqrt(a) / sqrt(2 pi), and the derivative
# multiplies it by i omega / samplingInterval (eq-marco-derivada,
# descomposicion.model.es.qmd:93-97), in rad/s. Both are returned on the
# record, samples by scales.
.transformCWT <- function(signal, samplingInterval, voiceCount, morletCenter) {
    n <- length(signal)
    Size <- 2^(1 + round(log2(n)))
    Left <- floor((Size - n) / 2)
    Right <- Size - n - Left
    Spectrum <- stats::fft(signal[c(Left:1, seq_len(n), n:(n - Right + 1))])
    AngularFrequency <- 2 * pi / Size * c(0:(Size / 2), (-Size / 2 + 1):(-1))
    Scales <- 2^(seq_len((log2(Size) - 1) * voiceCount) / voiceCount)
    Filters <- .evaluateMorlet(frequency = outer(AngularFrequency, Scales),
                               morletCenter = morletCenter) *
        rep(sqrt(Scales) / sqrt(2 * pi), each = Size)
    Observed <- Left + seq_len(n)
    list(coefficients = stats::mvfft(Filters * Spectrum, inverse = TRUE)[
             Observed, , drop = FALSE] / Size,
         derivative = stats::mvfft(Filters * (1i * AngularFrequency /
                                                  samplingInterval * Spectrum),
                                   inverse = TRUE)[Observed, , drop = FALSE] / Size)
}
