# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# DTW distance between the signal and each mode (eq-seleccion-dtw-local,
# eq-seleccion-dtw-recurrencia, eq-seleccion-dtw): the square root of the
# minimal accumulated cost G[N, N], without normalization by the path length.
# The recurrence G[i, j] = d[i, j] + min(G[i - 1, j - 1], G[i - 1, j],
# G[i, j - 1]), with G[0, 0] = 0 and infinite cost outside the matrix, is
# evaluated by anti-diagonals i + j = s, whose cells depend only on the two
# previous anti-diagonals; position i + 1 of each one holds row i.
measureDTW <- function(signal, modes) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || length(signal) < 1L ||
        !all(is.finite(signal)) || !is.matrix(modes) || !is.numeric(modes) ||
        !all(is.finite(modes)) || nrow(modes) != length(signal)) {
        stop("'signal' must be a nonempty finite numeric vector and 'modes' a ",
             "finite numeric matrix with one row per sample")
    }
    signal <- as.vector(signal)
    attributes(modes) <- list(dim = dim(modes))
    n <- length(signal)
    OUT <- numeric(ncol(modes))
    for (k in seq_len(ncol(modes))) {
        b <- modes[, k]
        Previous <- c(0, rep(Inf, n))
        Current <- rep(Inf, n + 1L)
        for (s in seq.int(2L, 2L * n)) {
            i <- seq.int(max(1L, s - n), min(n, s - 1L))
            AUX <- rep(Inf, n + 1L)
            AUX[i + 1L] <- (signal[i] - b[s - i])^2 +
                pmin(Previous[i], Current[i], Current[i + 1L])
            Previous <- Current
            Current <- AUX
        }
        OUT[k] <- sqrt(Current[n + 1L])
    }
    OUT
}
