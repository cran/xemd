# Directions of the vector sifting (memd.model:110-136): the Hammersley points
# i = 1, ..., D with the radical inverses rho_b(i) = sum_l a_l b^(-l-1) of the
# digits of i (memd.model:122-132) in the bases of the first C - 1 primes, in
# reverse order, and the coordinate i / D last, mapped by 2 p - 1 to b in the
# cube [-1, 1]^C. Each point becomes a unit vector through the angles of
# get.projections of emdr 1.2.3 (Pierre Masselot, revision
# a802208b98e506071bb4f3c7eea4d2e5ad0a16f7, GPL-3; memd.implementation:15-36),
# evaluated on that point alone: theta_j = atan2(|b_1..b_{C-j}|, b_{C+1-j}) and
# v = (cos theta_1, sin theta_1 cos theta_2, ..., sin theta_1 ... sin theta_{C-1}).
# Those D vectors cover a hemisphere; rows D + 1, ..., 2D are their antipodes in
# the same order, so that the maxima of the projections on rows d and d + D are
# the maxima and the minima of one projection.
.buildDirections <- function(directionCount, channelCount) {
    Primes <- integer()
    k <- 2L
    while (length(Primes) < channelCount - 1L) {
        if (all(k %% Primes != 0L)) Primes <- c(Primes, k)
        k <- k + 1L
    }
    IDX <- seq_len(directionCount)
    Points <- matrix(vapply(rev(Primes), function(base) {
        OUT <- numeric(directionCount)
        i <- IDX
        Scale <- 1 / base
        while (any(i > 0L)) {
            OUT <- OUT + i %% base * Scale
            i <- i %/% base
            Scale <- Scale / base
        }
        OUT
    }, numeric(directionCount)), nrow = directionCount)
    Points <- 2 * cbind(Points, IDX / directionCount) - 1
    Directions <- t(apply(Points, 1L, function(b) {
        C <- length(b)
        Theta <- atan2(sqrt(cumsum(b[-C]^2))[(C - 1L):1L], b[C:2L])
        Sine <- cumprod(sin(Theta))
        c(cos(Theta) * c(1, Sine[-(C - 1L)]), Sine[C - 1L])
    }))
    rbind(Directions, -Directions)
}
