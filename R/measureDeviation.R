# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Deviation sigma_h = |m / a| of the envelope mean m relative to the envelope
# amplitude a at every sample (emd.model:169-181); the caller forms m and a.
# Where a = 0 the ratio is undefined when m = 0, returned as NA, and +Inf when
# m != 0, as in emd.m and emd 0.8.1 (emd.rilling.implementation:189;
# sift.implementation:5).
.measureDeviation <- function(mean, amplitude) {
    OUT <- abs(mean / amplitude)
    OUT[mean == 0 & amplitude == 0] <- NA_real_
    OUT
}
