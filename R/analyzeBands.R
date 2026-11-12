# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Analysis coefficients of EWT, eq-ewt-analisis (ewt.sst.model.es.qmd:97-101):
# W_j = IDFT(conj(H_j) X) on the extension for each column of the bank
# (descomposicion.model.es.qmd:123-127). The bank is real and even in omega
# and the extension real, so W_j is real up to rounding.
.analyzeBands <- function(spectrum, filters) {
    Re(stats::mvfft(filters * spectrum, inverse = TRUE)) / nrow(filters)
}
