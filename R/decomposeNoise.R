# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Modes eta_k^(m) = E_k(w^(m)) of each noise realization (ceemdan.model:29,
# :260), decomposed with the configuration of the signal: a list of M
# matrices samples x K_m, in the order of the columns of noise.
.decomposeNoise <- function(noise, sifting) {
    lapply(seq_len(ncol(noise)), function(m) decomposeEMD(noise[, m], sifting = sifting)$modes)
}
