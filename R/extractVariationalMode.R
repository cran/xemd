# Extraction of one mode of SVMD (vmd.extensiones.es.qmd:316-350) from the
# spectrum of x on the grid of |omega| in rad/s, with the centers of the modes
# already extracted fixed and an initial center, both in rad/s. The penalty
# starts at alpha and doubles (eq-svmd-continuacion) while below maxAlpha; when
# doubling would exceed maxAlpha, one last computation runs at maxAlpha, and
# alpha = maxAlpha gives a single computation. Each penalty restarts the inner
# iteration with the spectrum and center reached, a zero multiplier and the
# penalty of the previous centers at the current alpha, and its iteration stops
# when the relative change (eq-svmd-parada-interna) falls below tolerance or
# after maxIterations. iterations adds the inner iterations of all penalties;
# converged requires every inner loop to meet the change test.
.extractVariationalMode <- function(spectrum, angularFrequency, previousCenters,
                                    initialCenter, alpha, maxAlpha, tau,
                                    tolerance, maxIterations) {
    State <- list(modeSpectrum = complex(length(spectrum)),
                  center = initialCenter)
    Iterations <- 0L
    Converged <- TRUE
    repeat {
        Penalty <- .computeCenterPenalty(angularFrequency = angularFrequency,
                                         angularCenters = previousCenters,
                                         alpha = alpha)
        State$multiplier <- complex(length(spectrum))
        for (Iteration in seq_len(maxIterations)) {
            Previous <- State$modeSpectrum
            State <- .stepSVMD(spectrum = spectrum, modeSpectrum = Previous,
                               center = State$center,
                               multiplier = State$multiplier,
                               angularFrequency = angularFrequency,
                               alpha = alpha, tau = tau, penalty = Penalty)
            # The mode is the case K = C = 1 of the change over modes and
            # channels.
            Change <- .measureSweepChange(current = matrix(State$modeSpectrum),
                                          previous = matrix(Previous))
            if (Change < tolerance) break
        }
        Iterations <- Iterations + Iteration
        Converged <- Converged && Change < tolerance
        if (alpha == maxAlpha) break
        alpha <- min(2 * alpha, maxAlpha)
    }
    list(modeSpectrum = State$modeSpectrum, center = State$center,
         iterations = Iterations, converged = Converged)
}
