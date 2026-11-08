# Time grid of a group of decomposeTSL: time, the N >= 2 increasing times in
# seconds that its OCID share, and series, the N x C matrix of their samples.
# The grid is uniform when diff(time) equals the mean step
# dt = (t_N - t_1) / (N - 1) under all.equal, the tolerance that gmsp applies
# to the rationality of the step (lib/R/AT2TS.R:201); a uniform grid returns
# unchanged with samplingInterval dt, rational or not. Otherwise dt becomes
# dt_r = .rationalizeStep(dt), the grid t_1 + k dt_r with
# k = 0, ..., floor((t_N - t_1) / dt_r + 1e-12), and each column the
# monotone cubic spline stats::splinefun(method = "monoH.FC") of (t - t_1, s)
# evaluated at k dt_r. Adapted from .regularize of gmsp 0.8.0 (commit
# 28e2ca2bf177), lib/R/regularize.R:7-39, Copyright (c) 2026 Alejandro Verri
# Kozlowski, MIT License, which starts from the smallest difference of the
# times and returns the grid from zero; here the step is the mean and the
# grid keeps the origin t_1.
.regularizeSeries <- function(time, series) {
    N <- length(time)
    Step <- (time[N] - time[1L]) / (N - 1)
    if (isTRUE(all.equal(diff(time), rep(Step, N - 1L)))) {
        return(list(time = time, series = series, samplingInterval = Step, regularized = FALSE))
    }
    Step <- .rationalizeStep(Step)
    Offset <- (0:floor((time[N] - time[1L]) / Step + 1e-12)) * Step
    list(time = time[1L] + Offset,
         series = apply(series, 2L, function(s) {
             stats::splinefun(time - time[1L], s, method = "monoH.FC")(Offset)
         }),
         samplingInterval = Step, regularized = TRUE)
}
