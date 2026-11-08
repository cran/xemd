# Adapted from LMD 1.2.1, revision ec11b721ec7b05f0137eeb67a619b28f62ef85a2.
# Author: Shubhra Prakash. Source license: Apache License (>= 2).
# Modified for xemd on 2026-10-11.

# Local mean and magnitude levels of LMD (eq-lmd-media-pareja,
# eq-lmd-amplitud-pareja; lmd.model:7-19). The nodes are the interior extrema
# of .findExtrema and the samples (1, x_1) and (N, x_N); consecutive nodes
# i, i + 1 give m_i = (n_i + n_{i+1}) / 2 and a_i = |n_i - n_{i+1}| / 2. A
# sample between two nodes takes the levels of their pair, a node on a sample
# the mean of the levels of its two pairs, and the first and last samples the
# levels of the first and last pairs, as local_mean_and_envelope of LMD 1.2.1.
# A plateau extremum may sit at a half-integer position, so a sample finds its
# pair by position. halfWidth is h = (K - 1) / 2 of the odd smoothing width
# K = max{3, floor(g / 3)}, plus one when even, with g the largest gap between
# consecutive nodes (lmd.model:49).
.buildLocalLevels <- function(signal, extrema) {
    N <- length(signal)
    Nodes <- rbind(c(1, signal[1L]), extrema$maxima, extrema$minima, c(N, signal[N]))
    Nodes <- Nodes[order(Nodes[, 1L]), , drop = FALSE]
    Position <- Nodes[, 1L]
    Value <- Nodes[, 2L]
    Last <- length(Position)
    Mean <- (Value[-Last] + Value[-1L]) / 2
    Amplitude <- abs(Value[-Last] - Value[-1L]) / 2
    Pair <- findInterval(seq_len(N), Position, rightmost.closed = TRUE)
    Levels <- list(mean = Mean[Pair], amplitude = Amplitude[Pair])
    IDX <- which(Position[-c(1L, Last)] %% 1 == 0) + 1L
    Levels$mean[Position[IDX]] <- (Mean[IDX - 1L] + Mean[IDX]) / 2
    Levels$amplitude[Position[IDX]] <- (Amplitude[IDX - 1L] + Amplitude[IDX]) / 2
    Width <- max(3, floor(max(diff(Position)) / 3))
    if (Width %% 2 == 0) Width <- Width + 1
    c(Levels, list(halfWidth = as.integer((Width - 1) / 2)))
}
