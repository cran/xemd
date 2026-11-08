# Rational sampling interval 1 / Fs, with Fs the smallest 2^a 5^b not below
# 1 / step once 1 / step is lowered by 16 machine epsilons relative to
# max(1, 1 / step); Fs is stored as an integer when it fits. Adapted from
# .rationalize of gmsp 0.8.0 (commit 28e2ca2bf177), lib/R/utils.R:125-140,
# Copyright (c) 2026 Alejandro Verri Kozlowski, MIT License: the same
# arithmetic, with the input guard left to the caller, whose steps are finite
# and positive.
.rationalizeStep <- function(step) {
    f <- 1 / step
    f <- f - 16 * .Machine$double.eps * max(1, f)
    b <- 0L:max(0L, ceiling(log(f, 5)) + 6L)
    a <- pmax(0L, ceiling(log(f, 2) - b * log(5, 2)))
    AUX <- 2^a * 5^b
    AUX <- AUX[AUX >= f]
    if (!length(AUX)) stop("no rate 2^a 5^b reaches 1 / step = ", 1 / step)
    Rate <- min(AUX)
    if (Rate <= .Machine$integer.max) Rate <- as.integer(round(Rate))
    1 / Rate
}
