# SPDX-FileCopyrightText: 2026 Alejandro Daniel Verri <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of xemd: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# Metrics by mode and by series of decomposition, the list of tables of
# decomposeTSL, for any of its methods. A mode c of a series x with N samples
# every dt = (t_N - t_1) / (N - 1) seconds, the step of the grid of its group,
# has: frequency, the power centroid of eq-vmd-centro on the reflection grid
# of 2N samples of the VMD family (.transformReflection, .computeCenter), in
# Hz, equal to the centers of VMD and SVMD but a center fixed at zero
# (vmd.model.es.qmd:271) and to the centroid of each channel, not the common
# center, in MVMD; bandwidth, the power deviation around it on the same grid,
# the square root of eq-vmd-banda-espectral over the energy of the mode
# (vmd.model.es.qmd:47-56), in Hz; energy, dt sum(c^2)
# (eq-marco-energia-discreta), in squared units of the series times seconds;
# energyFraction, sum(c^2) / sum(x^2), the E[k]/Eo of gmsp as a ratio (gmsp
# TS2IMF.R:273-308); and meanPeriod, N dt over the maxima of .findExtrema, or
# N dt without maxima, in seconds (modos.validacion.es.qmd:115, as
# assessWhiteNoise). A mode of zero energy determines no centroid
# (vmd.model.es.qmd:209): its frequency and bandwidth are NA. A series has
# modeCount; residueFraction, sum(r^2) / sum(x^2) with r its residue or
# fitResidual, zero when the method returns neither (EWT); and crossFraction,
# 1 - (sum_k sum(c_k^2) + sum(r^2)) / sum(x^2), the cross products of
# eq-marco-balance-energia over the energy of the series. A series of zero
# energy leaves its fractions undefined and signals xemdZeroEnergy. Returns
# the data.tables metrics (keys, mode, frequency, bandwidth, energy,
# energyFraction, meanPeriod) and metricsSeries (keys, modeCount,
# residueFraction, crossFraction), whose keys are those of the decomposition,
# its by columns, metadata, OCID and ID, in the order of its modes and of its
# series. decomposition is not modified.
measureModesTSL <- function(decomposition) {
    Keys <- .checkDecomposition(decomposition = decomposition,
                                columns = c("frequency", "bandwidth", "energy", "energyFraction", "meanPeriod",
                                            "modeCount", "residueFraction", "crossFraction"))
    Series <- unique(decomposition$signal[, Keys, with = FALSE])
    # Sums of the squared samples of a table by series, in the order of Series.
    sumSquares <- function(table) {
        IDX <- factor(Series[table, on = Keys, which = TRUE], levels = seq_len(nrow(Series)))
        vapply(split(table$s^2, IDX), sum, numeric(1L))
    }
    Energy <- sumSquares(decomposition$signal)
    if (any(Energy == 0)) {
        Label <- paste(Keys, vapply(Series[which.max(Energy == 0)], as.character, ""), sep = " = ", collapse = ", ")
        stop(errorCondition("the series has zero energy, so its energy fractions are undefined", unit = Label,
                            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    Modes <- unique(decomposition$modes[, c(Keys, "mode"), with = FALSE])
    IDX <- factor(Modes[decomposition$modes, on = c(Keys, "mode"), which = TRUE], levels = seq_len(nrow(Modes)))
    Times <- split(decomposition$modes$t, IDX)
    Samples <- split(decomposition$modes$s, IDX)
    Values <- vapply(seq_len(nrow(Modes)), function(k) {
        Mode <- Samples[[k]][order(Times[[k]])]
        n <- length(Mode)
        Step <- (max(Times[[k]]) - min(Times[[k]])) / (n - 1)
        Squares <- sum(Mode^2)
        Spectral <- c(NA_real_, NA_real_)
        if (Squares > 0) {
            Transform <- .transformReflection(signal = Mode, samplingInterval = Step)
            Power <- Mod(Transform$spectrum[, 1L])^2
            Center <- .computeCenter(power = Power, angularFrequency = Transform$angularFrequency)
            Spectral <- c(Center, sqrt(sum((Transform$angularFrequency - Center)^2 * Power) / sum(Power))) / (2 * pi)
        }
        c(Spectral, Step * Squares, Squares, n * Step / max(nrow(.findExtrema(Mode)$maxima), 1L))
    }, numeric(5L))
    Owner <- Series[Modes, on = Keys, which = TRUE]
    Remainder <- intersect(c("residue", "fitResidual"), names(decomposition))
    Rest <- if (length(Remainder)) sumSquares(decomposition[[Remainder]]) else numeric(nrow(Series))
    list(metrics = cbind(Modes, data.table(frequency = Values[1L, ], bandwidth = Values[2L, ],
                                           energy = Values[3L, ], energyFraction = Values[4L, ] / Energy[Owner],
                                           meanPeriod = Values[5L, ])),
         metricsSeries = cbind(Series, data.table(modeCount = tabulate(Owner, nbins = nrow(Series)),
                                                  residueFraction = Rest / Energy,
                                                  crossFraction = 1 - (sumSquares(decomposition$modes) + Rest) /
                                                      Energy)))
}
