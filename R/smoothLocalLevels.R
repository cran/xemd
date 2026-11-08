# Repeated uniform moving average of the local levels of LMD
# (eq-lmd-promedio-discreto, eq-lmd-semiancho-borde, eq-lmd-promedio-borde;
# lmd.model:25-47), with the half width h of the levels in place of the 2 of
# the five-point example: one pass replaces sample n by the mean of the
# 2 h_n + 1 samples around it, h_n = min{h, n - 1, N - n}, so the first and
# last samples stay fixed. Passes repeat while two consecutive samples of the
# result are equal, compared exactly as the source states no tolerance
# (lmd.model:47), up to maxSmoothPasses passes; the last pass is returned.
.smoothLocalLevels <- function(levels, halfWidth, maxSmoothPasses) {
    N <- length(levels)
    Reach <- pmin(halfWidth, seq_len(N) - 1L, N - seq_len(N))
    for (Pass in seq_len(maxSmoothPasses)) {
        levels <- vapply(seq_len(N), function(n) {
            sum(levels[(n - Reach[n]):(n + Reach[n])]) / (2 * Reach[n] + 1)
        }, numeric(1L))
        if (!any(diff(levels) == 0)) break
    }
    levels
}
