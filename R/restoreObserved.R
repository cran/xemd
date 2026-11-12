# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Inverse of .transformReflection for each column: inverse DFT divided by T
# (eq-marco-idft) at the observed positions of the extension
# (vmd.model.es.qmd:269). The spectra are Hermitian, so the discarded
# imaginary part is rounding.
.restoreObserved <- function(spectrum, observed) {
    Re(stats::mvfft(spectrum, inverse = TRUE)[observed, , drop = FALSE]) /
        nrow(spectrum)
}
