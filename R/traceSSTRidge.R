# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Adapted from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979:
# curve_ext of the Synchrosqueezing Toolbox, original author Jianfeng Lu,
# maintainer Eugene Brevdo; Fortran port: Dongik Jang.
# SPDX-License-Identifier: LGPL-2.0-or-later
# xemd distributes this file under GPL-3, as section 3 of the GNU Library
# General Public License version 2 permits.
# Modified for xemd on 2026-10-15.

# Ridge of a synchrosqueezed transform, samples by frequency cells, by the
# dynamic programme of curve_ext in the Synchrosqueezing Toolbox (supplement
# sst, eq-deferred-sst-energy and eq-deferred-sst-dynamic): the cost of a cell
# is E = -log(|T|^2 / Q + 1e-8), Q the energy of the whole transform; with
# z = log2(frequency), F_{j,1} = E_{j,1} and F_{j,n} = E_{j,n} +
# min_k [F_{k,n-1} + penalty (z_j - z_k)^2] over every cell k. The last index
# is the first minimizer of F_{.,N}, and each predecessor is the first cell
# whose residual in the recurrence is below 1e-8 in modulus. Returns one cell
# index per sample. A transform of zero energy leaves Q undefined; a sample
# without predecessor leaves the trajectory undefined.
.traceSSTRidge <- function(transform, frequency, penalty) {
    Power <- Re(transform)^2 + Im(transform)^2
    if (sum(Power) == 0) {
        stop(errorCondition("the transform has zero energy", class = c("xemdZeroEnergy", "moduleFailure")))
    }
    Cost <- -log(Power / sum(Power) + 1e-8)
    Jump <- penalty * outer(log2(frequency), log2(frequency), "-")^2
    Accumulated <- matrix(0, nrow = ncol(Cost), ncol = nrow(Cost))
    Accumulated[, 1L] <- Cost[1L, ]
    for (m in seq_len(nrow(Cost))[-1L]) {
        Accumulated[, m] <- Cost[m, ] + apply(Accumulated[, m - 1L] + Jump, 2L, min)
    }
    Path <- integer(nrow(Cost))
    Path[nrow(Cost)] <- which.min(Accumulated[, nrow(Cost)])
    for (m in rev(seq_len(nrow(Cost) - 1L))) {
        Next <- Path[m + 1L]
        Residual <- Accumulated[Next, m + 1L] - Cost[m + 1L, Next] -
            Accumulated[, m] - Jump[, Next]
        Path[m] <- which(abs(Residual) < 1e-8)[1L]
        if (is.na(Path[m])) {
            stop("no cell at sample ", m, " precedes cell ", Next, " of the ridge")
        }
    }
    Path
}
