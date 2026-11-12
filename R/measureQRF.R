# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Quality of reconstruction factor of an estimate against its reference
# (eq-validacion-qrf), in dB; Inf when both coincide.
measureQRF <- function(estimate, reference) {
    if (!is.numeric(estimate) || !is.null(dim(estimate)) || !all(is.finite(estimate)) ||
        !is.numeric(reference) || !is.null(dim(reference)) ||
        !all(is.finite(reference)) || length(estimate) != length(reference)) {
        stop("'estimate' and 'reference' must be finite numeric vectors of equal length")
    }
    estimate <- as.vector(estimate)
    reference <- as.vector(reference)
    Norm <- sqrt(sum(reference^2))
    if (Norm == 0) stop("'reference' must have a nonzero norm")
    20 * log10(Norm / sqrt(sum((reference - estimate)^2)))
}
