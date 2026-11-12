# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Fractional Gaussian noise diagnostic of the IMF of a record of 2048 samples
# with Hurst exponent hurst, at the tabulated level probability (supplement
# validation ruido-fraccional, eq-validacion-ruido-energia-inicial,
# -calibrado, -razon and -umbral-explicito; table 3.1 of Flandrin, Goncalves
# and Rilling). The energy of IMF k is W[k] = sum c_k^2; the first IMF gives
# the noise scale, and from the second on the reference is
# W[1] / beta_H rho_H^(-2 (1 - H) k), rho_H = 2.01 + 0.2 (H - 0.5) +
# 0.12 (H - 0.5)^2, and the threshold W_H[k] 2^(2^(a k + b)). An IMF lies
# outside when its energy exceeds the threshold. The coefficients exist only
# for the exponents, levels and length of the table, and are not
# interpolated. The first IMF is not evaluated: NA in reference, upper and
# outside. A first IMF of zero energy leaves the reference undefined.
assessFractionalNoise <- function(modes, hurst, probability) {
    HURST <- c(0.2, 0.5, 0.8)
    LEVEL <- c(0.95, 0.99)
    BETA <- c(0.487, 0.719, 1.025)
    SLOPE <- rbind(c(0.458, 0.452), c(0.474, 0.460), c(0.497, 0.495))
    INTERCEPT <- rbind(c(-2.435, -1.951), c(-2.449, -1.919), c(-2.331, -1.833))
    stopifnot(
        "modes must be a finite numeric matrix of 2048 samples by at least two IMF" =
            is.numeric(modes) && is.matrix(modes) && nrow(modes) == 2048L &&
            ncol(modes) >= 2L && all(is.finite(modes)),
        "hurst must be 0.2, 0.5 or 0.8" =
            is.numeric(hurst) && length(hurst) == 1L && hurst %in% HURST,
        "probability must be 0.95 or 0.99" =
            is.numeric(probability) && length(probability) == 1L && probability %in% LEVEL)
    attributes(modes) <- list(dim = dim(modes))
    i <- match(hurst, HURST)
    j <- match(probability, LEVEL)
    Energy <- colSums(modes^2)
    if (Energy[1L] == 0) {
        stop(errorCondition("the first IMF has zero energy", class = c("xemdZeroEnergy", "moduleFailure")))
    }
    k <- seq_along(Energy)
    Ratio <- 2.01 + 0.2 * (hurst - 0.5) + 0.12 * (hurst - 0.5)^2
    Reference <- c(NA, (Energy[1L] / BETA[i] * Ratio^(-2 * (1 - hurst) * k))[-1L])
    Upper <- Reference * 2^(2^(SLOPE[i, j] * k + INTERCEPT[i, j]))
    list(energy = Energy, reference = Reference, upper = Upper,
         outside = Energy > Upper)
}
