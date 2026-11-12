# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Penalty of SVMD on the modes already extracted
# (eq-svmd-penalizacion-espectral, vmd.extensiones.es.qmd:237-244): at each
# |omega| of the grid, the sum over the previous centers of
# 1 / (alpha^2 (omega - omega_i)^4), with the current alpha in s^2 and the
# centers in rad/s. It is zero without previous centers and +Inf at a bin equal
# to a previous center, where the update of the mode vanishes
# (vmd.extensiones.es.qmd:253).
.computeCenterPenalty <- function(angularFrequency, angularCenters, alpha) {
    rowSums(1 / (alpha * outer(angularFrequency, angularCenters, "-")^2)^2)
}
