# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979.
# MATLAB original: Eugene Brevdo; R port: Dongik Jang, Hee-Seok Oh and Donghoh Kim.
# Synchrosqueezing toolbox attribution also names Hau-Tieng Wu.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-06 and 2026-10-11.

# Fourier transform of the Morlet wavelet with the mean correction, so that
# psi^(0) = 0, at frequencies in rad/sample times the scale
# (synchrosqueezing.implementation.es.qmd:3): with mu = morletCenter,
# c pi^(-1/4) (exp(-(mu - xi)^2 / 2) - exp(-mu^2 / 2) exp(-xi^2 / 2)), where
# c = (1 + exp(-mu^2) - 2 exp(-3 mu^2 / 4))^(-1/2) makes the integral of its
# square over the real line one.
.evaluateMorlet <- function(frequency, morletCenter) {
    Normalization <- (1 + exp(-morletCenter^2) -
                          2 * exp(-3 / 4 * morletCenter^2))^(-1 / 2)
    Normalization * pi^(-1 / 4) * (exp(-(morletCenter - frequency)^2 / 2) -
        exp(-morletCenter^2 / 2) * exp(-frequency^2 / 2))
}
