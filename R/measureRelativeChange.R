# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Relative change between consecutive iterates of iterative filtering,
# SD_q = ||f^(q) - f^(q-1)||_2 / ||f^(q-1)||_2 (eq-if-parada-relativa;
# lmd.model:290-297): a ratio of norms, not of squared norms. By Parseval it
# takes the same value on two real states or on their full discrete spectra.
# A null previous state leaves it undefined (NaN).
.measureRelativeChange <- function(previous, current) {
    sqrt(sum(Mod(current - previous)^2)) / sqrt(sum(Mod(previous)^2))
}
