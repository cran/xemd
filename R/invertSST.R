# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979.
# MATLAB original: Eugene Brevdo; R port: Dongik Jang, Hee-Seok Oh and Donghoh Kim.
# Synchrosqueezing toolbox attribution also names Hau-Tieng Wu.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-11.

# Inversion of the synchrosqueezed transform, eq-sst-inversion-real
# (ewt.sst.model.es.qmd:177-181), discretized as SynchWave
# (synchrosqueezing.implementation.es.qmd:95): x(b) = Re sum_f T(b, f) / C
# with C = C_psi / (sqrt(2 pi) 2 ln 2). C_psi is eq-sst-constante
# (ewt.sst.model.es.qmd:161-166) for the real Morlet wavelet, integrated by
# pracma::quadgk after the change xi = (x / (1 - x))^2 with SynchWave's
# tolerance sqrt(eps); the other factors are the 2 of the real part, the ln 2
# of the logarithmic scale measure ln 2 / voiceCount, whose 1 / voiceCount T
# already carries, and the 1 / sqrt(2 pi) of the CWT filter. The recovery is
# approximate.
.invertSST <- function(transform, morletCenter) {
    integrand <- function(x) {
        Root <- x / (1 - x)
        2 * .evaluateMorlet(frequency = Root^2, morletCenter = morletCenter) /
            (Root * (1 - x)^2)
    }
    Admissibility <- pracma::quadgk(integrand, a = 0, b = 1,
                                    tol = sqrt(.Machine$double.eps))
    Re(rowSums(transform)) * sqrt(2 * pi) * 2 * log(2) / Admissibility
}
