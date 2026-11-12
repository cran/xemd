# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Fault correlation factor of the signal with each mode (eq-seleccion-fcf): a
# normalized inner product without removal of means.
measureFCF <- function(signal, modes) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal)) ||
        !is.matrix(modes) || !is.numeric(modes) || !all(is.finite(modes)) ||
        nrow(modes) != length(signal)) {
        stop("'signal' must be a finite numeric vector and 'modes' a finite ",
             "numeric matrix with one row per sample")
    }
    signal <- as.vector(signal)
    attributes(modes) <- list(dim = dim(modes))
    Energy.signal <- sum(signal^2)
    Energy.modes <- colSums(modes^2)
    if (Energy.signal == 0 || any(Energy.modes == 0)) {
        stop(errorCondition("'signal' and every column of 'modes' must have a nonzero norm",
                            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    as.vector(crossprod(modes, signal)) / sqrt(Energy.signal * Energy.modes)
}
