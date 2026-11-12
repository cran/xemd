# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Spectral grid of the VMD family. The record of N samples is extended to
# T = 2N samples by reflecting ceiling(N / 2) samples at the start and the
# remaining ones at the end (vmd.model.es.qmd:269), and every column is
# transformed with the unnormalized DFT (eq-marco-dft). The angular grid is
# |omega| in rad/s, 2 pi min(j, T - j) / (T dt) for j = 0, ..., T - 1 in the
# order of stats::fft (descomposicion.model.es.qmd:63, 135). observed holds
# the positions of the record inside the extension.
.transformReflection <- function(signal, samplingInterval) {
    signal <- matrix(signal, nrow = NROW(signal))
    n <- nrow(signal)
    Left <- ceiling(n / 2)
    Positions <- c(Left + 1 - seq_len(Left), seq_len(n),
                   n + 1 - seq_len(n - Left))
    j <- seq_len(2L * n) - 1L
    list(spectrum = stats::mvfft(signal[Positions, , drop = FALSE]),
         angularFrequency = 2 * pi * pmin(j, 2L * n - j) /
             (2 * n * samplingInterval),
         observed = Left + seq_len(n))
}
