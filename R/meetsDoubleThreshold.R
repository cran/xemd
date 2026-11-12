# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Double threshold on sigma_h (emd.model:183-212): the fraction of the N
# samples with sigma_h >= theta_1 is at most alpha and no sample has
# sigma_h >= theta_2. Undefined samples (NA) are no exceedance and stay in the
# denominator N. The count condition |N_e - N_z| <= 1 belongs to the caller.
.meetsDoubleThreshold <- function(deviation, thresholds, fraction) {
    sum(deviation >= thresholds[1L], na.rm = TRUE) / length(deviation) <= fraction &&
        !any(deviation >= thresholds[2L], na.rm = TRUE)
}
