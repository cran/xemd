# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Criterion SD of the 1998 formulation (emd.model:106-124): the sum over all N
# samples of (h_{j-1} - h_j)^2 / h_{j-1}^2, without the factor 1/N. A null
# denominator leaves an infinite or undefined term; it is not regularized.
.measureSD <- function(previous, current) {
    sum((previous - current)^2 / previous^2)
}
