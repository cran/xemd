# First vector mode E_1 by sifting (memd.model:54-62, :102-108): h_0 is the
# signal and h_{s+1} = h_s - m[h_s]. With envelope "trajectory" the mean is
# (1 / D) sum_d e_d over the envelopes of the trajectory at the maxima of its
# projections on the D rows of directions (memd.model:26-36, :84-100); with
# "projection" it is the mean of the projection scheme (.computeProjectedMean).
# Each state is examined once, in this order: the rule accepts it where the
# rule is defined (converged TRUE); s equals fixedSift under "fixed"
# (converged NA); s equals maxSift under an adaptive rule, or a direction has
# no envelope (converged FALSE); otherwise one subtraction. The adaptive rules
# act on the projections (Rehman and Mandic 2011, section II): "consecutive"
# accepts when every projection keeps |E - Z| <= 1 in consecutiveCount
# successive states after a subtraction, the run restarting when one does not;
# "doubleThreshold" accepts a state, the signal included, whose projections
# all keep |E - Z| <= 1 and whose ratio |m| / a, with the amplitude a of
# amplitudeRule, meets the double threshold of the scalar sifting. Returns
# a NULL mode when the initial mean does not exist. Failed attempts carry
# their cause and the first direction without an envelope. sifting is complete and
# validated by the public entry that receives it.
.extractVectorIMF <- function(signal, directions, envelope, sifting) {
    Rule <- sifting$stopRule
    Current <- signal
    Run <- 0L
    s <- 0L
    repeat {
        if (envelope == "projection") {
            Projected <- .computeProjectedMean(Current, directions = directions, boundary = sifting$boundary)
            Mean <- Projected$mean
            Direction <- Projected$direction
        }
        if (envelope == "trajectory") {
            Projections <- Current %*% t(directions)
            Projections <- split(Projections, col(Projections))
            Extrema <- lapply(Projections, .findExtrema)
            Envelopes <- Map(function(projection, extrema) {
                .buildVectorEnvelope(Current, projection = projection, extrema = extrema,
                                     boundary = sifting$boundary)
            }, Projections, Extrema)
            Mean <- if (!any(vapply(Envelopes, is.null, logical(1L)))) Reduce(`+`, Envelopes) / length(Envelopes)
            Direction <- unname(which(vapply(Envelopes, is.null, logical(1L)))[1L])
            Compatible <- all(vapply(Extrema, function(x) {
                abs(nrow(x$maxima) + nrow(x$minima) - x$crossings) <= 1L
            }, logical(1L)))
        }
        if (s == 0L && is.null(Mean)) {
            return(list(mode = NULL, iterations = s, converged = FALSE,
                        stopReason = "insufficientExtrema", direction = Direction))
        }
        Accepted <- FALSE
        if (Rule == "consecutive" && s >= 1L) {
            Run <- if (Compatible) Run + 1L else 0L
            Accepted <- Run >= sifting$consecutiveCount
        }
        if (Rule == "doubleThreshold" && !is.null(Mean)) {
            Amplitude <- .measureVectorAmplitude(Envelopes, mean = Mean, amplitudeRule = sifting$amplitudeRule)
            Deviation <- .measureDeviation(mean = sqrt(rowSums(Mean^2)), amplitude = Amplitude)
            Accepted <- Compatible &&
                .meetsDoubleThreshold(Deviation, thresholds = sifting$thresholds, fraction = sifting$fraction)
        }
        if (Accepted) return(list(mode = Current, iterations = s, converged = TRUE))
        if (Rule == "fixed" && s == sifting$fixedSift) {
            return(list(mode = Current, iterations = s, converged = NA))
        }
        if ((Rule != "fixed" && s == sifting$maxSift) || is.null(Mean)) {
            return(list(mode = Current, iterations = s, converged = FALSE,
                        stopReason = if (is.null(Mean)) "insufficientExtrema" else "iterationLimit",
                        direction = Direction))
        }
        Current <- Current - Mean
        s <- s + 1L
    }
}
