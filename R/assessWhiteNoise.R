# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# White-noise significance test of Wu and Huang on each column
# (modos.validacion:107-250). The energy E_k = sum(c_k^2) / N
# (eq-validacion-wu-energia) is compared with the bounds E_p of the
# probabilities p- < p+ at the mean period T_k, N over the count of maxima of
# .findExtrema, or N without maxima (:115). "gaussian" is the line of
# eq-validacion-wu-cuantil-log, ln E_p = -ln T + z_p sqrt(2 / N) T^(1/2);
# "chiSquare" is the quantile of eq-validacion-wu-cuantil-energia,
# E_p = q(p; r) / N with r = N / T (eq-validacion-wu-grados on the reference
# eq-validacion-wu-referencia). The bounds hold for components of a
# normalized signal, which the caller supplies (:145-151). With rescalingMode
# r the energy of component r is placed on its upper bound and every other
# energy is multiplied by E_p+(T_r) / E_r (:246-248). outside is
# E < lower or E > upper, both tails (:244).
assessWhiteNoise <- function(modes, probabilities, bound, rescalingMode) {
    if (!is.matrix(modes) || !is.numeric(modes) || !all(is.finite(modes))) {
        stop("'modes' must be a finite numeric matrix")
    }
    if (!is.numeric(probabilities) || length(probabilities) != 2L ||
        !all(is.finite(probabilities)) || probabilities[1L] <= 0 ||
        probabilities[1L] >= probabilities[2L] || probabilities[2L] >= 1) {
        stop("'probabilities' must be two numbers with 0 < probabilities[1] < probabilities[2] < 1")
    }
    if (!identical(bound, "gaussian") && !identical(bound, "chiSquare")) {
        stop("'bound' must be \"gaussian\" or \"chiSquare\"")
    }
    if (!is.null(rescalingMode) &&
        !(is.numeric(rescalingMode) && length(rescalingMode) == 1L &&
          rescalingMode %in% seq_len(ncol(modes)))) {
        stop("'rescalingMode' must be NULL or the index of a column of 'modes'")
    }
    attributes(modes) <- list(dim = dim(modes))
    n <- nrow(modes)
    Energy <- colSums(modes^2)
    if (any(Energy == 0)) {
        stop(errorCondition("'modes' has a column of zero energy", class = c("xemdZeroEnergy", "moduleFailure")))
    }
    Energy <- Energy / n
    Period <- vapply(seq_len(ncol(modes)), function(k) {
        Count <- nrow(.findExtrema(modes[, k])$maxima)
        if (Count == 0L) return(n)
        n / Count
    }, numeric(1L))
    computeBound <- function(p) {
        if (bound == "gaussian") return(exp(-log(Period) + stats::qnorm(p) * sqrt(2 / n) * sqrt(Period)))
        if (bound == "chiSquare") return(stats::qchisq(p, df = n / Period) / n)
    }
    Lower <- computeBound(probabilities[1L])
    Upper <- computeBound(probabilities[2L])
    if (!is.null(rescalingMode)) {
        Energy <- Energy * Upper[rescalingMode] / Energy[rescalingMode]
        Energy[rescalingMode] <- Upper[rescalingMode]
    }
    list(energy = Energy, period = Period, lower = Lower, upper = Upper,
         outside = Energy < Lower | Energy > Upper)
}
