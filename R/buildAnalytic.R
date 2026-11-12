# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Analytic signal z = c + i Hc of each column (eq-hilbert-analitica): the DFT
# of the column is multiplied by 1 + sgn of the bin frequency
# (eq-hilbert-unilateral). The zero bin and, for an even number of samples, the
# Nyquist bin keep the factor 1, which makes Re z equal the column. The DFT
# implies the periodic extension of the record.
.buildAnalytic <- function(signal) {
    n <- nrow(signal)
    k <- seq_len(n) - 1L
    Factor <- 1 + sign(n - 2 * k) * (k > 0L)
    stats::mvfft(Factor * stats::mvfft(signal), inverse = TRUE) / n
}
