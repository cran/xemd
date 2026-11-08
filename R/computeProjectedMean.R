# Mean of the projection scheme of bivariate EMD (memd.model:40-52): for each
# direction v_d = (cos phi_d, sin phi_d) the envelope of the maxima of the
# projection p_d = x v_d (memd.model:17-22), oriented by exp(i phi_d), which in
# the real coordinates of the plane is that scalar envelope times v_d; the sum
# over the D directions is scaled by 2 / D. Returns the mean and the first
# direction without an envelope, if any; that case has a NULL mean.
.computeProjectedMean <- function(signal, directions, boundary) {
    Mean <- 0
    for (d in seq_len(nrow(directions))) {
        Direction <- directions[d, ]
        Projection <- drop(signal %*% Direction)
        Envelopes <- .buildEnvelopes(Projection, extrema = .findExtrema(Projection), boundary = boundary)
        if (is.null(Envelopes)) return(list(mean = NULL, direction = d))
        Mean <- Mean + outer(Envelopes$upper, Direction)
    }
    list(mean = 2 * Mean / nrow(directions), direction = NA_integer_)
}
