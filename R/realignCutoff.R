# Adapted from pytvfemd 0.3.14 (Stefano Bianchi and Alessandro Longo), GPL-3;
# modified for xemd on 2026-10-11.
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2022-  Stefano Bianchi, Alessandro Longo

# Realigned cut-off of TVF-EMD in rad per sample (tvf.model:281-295). The
# cutoff above 0.9 pi (0.45 cycles per sample) is set to 0.9 pi and the samples
# where it is undefined or not positive are replaced by the pchip interpolation
# of the positive ones, the first and the last sample taking the first and the
# last positive value, as pytvfemd 0.3.14 treats the domain of the cutoff
# (tvf.model:279; pytvfemd.py:97-99, :213-232, :249-250). Then algorithm 2:
# with u_i the maxima of signal, an event e_j is
# declared at u_i when the cutoff over [u_i, u_{i+1}) has a positive minimum
# and (max - min) / min > threshold (eq-tvf-intermitencia); it rises or falls
# with the sign of cutoff(u_{i+1}) - cutoff(u_i), a half-integer position
# taking the mean of its two central samples. Between three neighbouring
# events a rising e_j makes [e_{j-1}, e_j) a low zone and a falling e_j makes
# [e_j, e_{j+1}) one; the low zones are replaced by the pchip interpolation of
# the rest (pytvfemd.py:231-232). Returns NULL when fewer than three samples
# have a positive cutoff, the nodes that pracma::pchip needs.
.realignCutoff <- function(cutoff, signal, threshold) {
    Sample <- seq_along(cutoff)
    cutoff <- pmin(cutoff, 0.9 * pi)
    Inside <- !is.na(cutoff) & cutoff > 0
    if (sum(Inside) < 3L) return(NULL)
    Border <- c(1L, length(cutoff))
    cutoff[Border] <- cutoff[range(which(Inside))]
    Inside[Border] <- TRUE
    if (!all(Inside)) cutoff[!Inside] <- pracma::pchip(Sample[Inside], cutoff[Inside], Sample[!Inside])
    Maxima <- .findExtrema(signal)$maxima[, "position"]
    Span <- findInterval(Sample, Maxima)
    MASK <- Span >= 1L & Span < length(Maxima)
    Lowest <- tapply(cutoff[MASK], Span[MASK], min)
    Highest <- tapply(cutoff[MASK], Span[MASK], max)
    IDX <- as.integer(names(Lowest))[Lowest > 0 & (Highest - Lowest) / Lowest > threshold]
    valueAt <- function(position) (cutoff[floor(position)] + cutoff[ceiling(position)]) / 2
    Events <- Maxima[IDX]
    Edges <- sign(valueAt(Maxima[IDX + 1L]) - valueAt(Events))
    Low <- logical(length(cutoff))
    for (j in seq_along(Events)[-c(1L, length(Events))]) {
        if (Edges[j] > 0) Low <- Low | (Sample >= Events[j - 1L] & Sample < Events[j])
        if (Edges[j] < 0) Low <- Low | (Sample >= Events[j] & Sample < Events[j + 1L])
    }
    if (any(Low)) cutoff[Low] <- pracma::pchip(Sample[!Low], cutoff[!Low], Sample[Low])
    cutoff
}
