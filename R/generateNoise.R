# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Perturbation realizations w^(m), m = 1, ..., realizationCount, of
# eq-eemd-blancura (eemd.model:27-33; ceemdan.model:29): independent Gaussian
# samples of zero mean and unit variance drawn by stats::rnorm from the
# caller's RNG stream, without a seed of their own. Column m is realization
# m, so the samples are drawn realization by realization.
.generateNoise <- function(sampleCount, realizationCount) {
    isCount <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
    if (!(isCount(sampleCount) && isCount(realizationCount))) {
        stop("sampleCount and realizationCount must be whole numbers of at least one")
    }
    matrix(stats::rnorm(sampleCount * realizationCount), nrow = sampleCount, ncol = realizationCount)
}
