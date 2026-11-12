# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Average by extraction order of aligned decompositions: D^0_j with the
# divisor M for "realizations" (eq-eemd-divisores, eemd.model:195-199) or
# D~_j with the divisor M_j of the realizations that have mode j for
# "available" (eq-eemd-media-disponibles, :189-193). The residue is the mean
# of the M residues with either divisor (eq-eemd-media-residuo, :57-61).
# iterations sums the subtractions of each mode (eq-eemd-trabajo, :299-305);
# converged combines the flags of the realizations that have the mode, TRUE
# only if all accepted, with FALSE over NA.
.averageDecompositions <- function(aligned, divisor) {
    K <- dim(aligned$modes)[2L]
    Divisor <- switch(divisor,
                      realizations = rep(length(aligned$modeCount), K),
                      available = vapply(seq_len(K), function(j) sum(aligned$modeCount >= j), integer(1L)),
                      stop("unrecognized divisor: ", divisor))
    list(modes = rowSums(aligned$modes, dims = 2L) / rep(Divisor, each = dim(aligned$modes)[1L]),
         residue = rowMeans(aligned$residue),
         iterations = as.integer(rowSums(aligned$iterations)),
         converged = vapply(seq_len(K), function(j) all(aligned$converged[j, aligned$modeCount >= j]),
                            logical(1L)))
}
