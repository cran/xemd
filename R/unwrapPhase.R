# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Continuous phase of each column of an analytic signal (eq-hilbert-polar).
# Within each run of samples with |z| > 0 the phase starts at the principal
# value Arg z and adds the increments Arg(z[n + 1] Conj(z[n])) taken in
# (-pi, pi], so an increment of exactly -pi counts as pi. A zero of z leaves
# its sample undefined (NA), and the next run starts again at its principal
# value.
.unwrapPhase <- function(analytic) {
    n <- nrow(analytic)
    Phase <- matrix(NA_real_, nrow = n, ncol = ncol(analytic))
    for (k in seq_len(ncol(analytic))) {
        z <- analytic[, k]
        Defined <- Mod(z) > 0
        Start <- Defined & c(TRUE, !Defined[-n])
        Step <- c(0, Arg(z[-1L] * Conj(z[-n])))
        Step[Step == -pi] <- pi
        Step[Start] <- Arg(z[Start])
        Phase[Defined, k] <- stats::ave(Step, cumsum(Start), FUN = cumsum)[Defined]
    }
    Phase
}
