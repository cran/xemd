# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# One sweep of the VMD recurrence (vmd.model.es.qmd:167-221) on C >= 1
# channels sharing K centers, as in MVMD (vmd.extensiones.es.qmd:48-141).
# spectrum and multiplier are T x C, modeSpectra T x K x C, angularCenters in
# rad/s. Mode k is updated after modes 1, ..., k - 1 and sees their new spectra
# (eq-vmd-entrada-iteracion, eq-mvmd-entrada); its denominator uses its center
# from the previous sweep and is common to the channels (eq-vmd-actualizacion,
# eq-mvmd-modo), and its new center weighs the power summed over channels
# (eq-vmd-centro, eq-mvmd-centro). The dual ascent of each channel closes the
# sweep (eq-vmd-dual, eq-mvmd-dual).
.sweepVMD <- function(spectrum, modeSpectra, angularCenters, multiplier,
                      angularFrequency, alpha, tau, fixedZeroCenter) {
    Total <- rowSums(aperm(modeSpectra, c(1L, 3L, 2L)), dims = 2L)
    for (k in seq_along(angularCenters)) {
        Others <- Total - modeSpectra[, k, ]
        modeSpectra[, k, ] <- (spectrum - Others + multiplier / 2) /
            (1 + 2 * alpha * (angularFrequency - angularCenters[k])^2)
        Total <- Others + modeSpectra[, k, ]
        # The fixed zero center is not updated (vmd.model.es.qmd:271), yet its
        # mode without energy is the same undefined case.
        Center <- .computeCenter(
            power = rowSums(Mod(modeSpectra[, k, , drop = FALSE])^2),
            angularFrequency = angularFrequency)
        if (!fixedZeroCenter || k > 1L) angularCenters[k] <- Center
    }
    list(modeSpectra = modeSpectra, angularCenters = angularCenters,
         multiplier = multiplier + tau * (spectrum - Total))
}
