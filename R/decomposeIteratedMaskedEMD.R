# Iterated masking EMD, itEMD (mascaras.model:182-229): a fixed point of the
# mask frequencies, f^(q+1) = F_x(f^(q)) (eq-mascaras-punto-fijo). Each
# iteration decomposes the signal with the uniform-phase masked decomposition
# of frequencies f^(q), phaseCount phases and amplitude epsilon_0 sd_N of the
# residue (eq-mascaras-coste-itemd), and updates each frequency to the mean
# instantaneous frequency of its mode weighted by the squared amplitude, over
# the samples where the frequency is defined (eq-mascaras-analitica,
# eq-mascaras-frecuencia-instantanea, eq-mascaras-actualizacion). The
# iteration stops with "tolerance" when every relative change
# |f^(q+1) - f^(q)| / f^(q) is below tolerance (eq-mascaras-cambio-relativo),
# or with "iterationLimit" after maxIterations updates. It returns the last
# decomposition with the frequencies, in Hz, that produced it, the count of
# updates and the termination of that decomposition. Defaults: iter_th = 0.1
# and max_iter = 15 of iterated_mask_sift in emd 0.8.1
# (mascaras.iterativas.implementation:107), the code of Andrew Quinn, coauthor
# of itEMD (mascaras.implementation:289).
decomposeIteratedMaskedEMD <- function(signal, samplingInterval, initialFrequencies,
                                       phaseCount, amplitudeRatio, tolerance = 0.1,
                                       maxIterations = 15, sifting) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal))) {
        stop("signal must be a numeric vector of finite real values without dim")
    }
    isPositive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
    if (!isPositive(samplingInterval)) {
        stop("samplingInterval must be a finite number of seconds greater than zero")
    }
    if (!(is.numeric(initialFrequencies) && is.null(dim(initialFrequencies)) &&
          length(initialFrequencies) >= 1L && all(is.finite(initialFrequencies)) &&
          all(initialFrequencies > 0) && all(initialFrequencies <= 1 / (2 * samplingInterval)))) {
        stop("initialFrequencies must be frequencies in Hz in (0, 1 / (2 samplingInterval)]")
    }
    if (!(isPositive(phaseCount) && phaseCount >= 2 && phaseCount == round(phaseCount))) {
        stop("phaseCount must be a whole number of at least two")
    }
    if (!isPositive(amplitudeRatio)) {
        stop("amplitudeRatio must be a finite number greater than zero")
    }
    if (!isPositive(tolerance)) stop("tolerance must be a finite number greater than zero")
    if (!(isPositive(maxIterations) && maxIterations == round(maxIterations))) {
        stop("maxIterations must be a whole number of at least one")
    }
    sifting <- .checkSifting(sifting)
    if (!sifting$stopRule %in% c("SD", "Snumber", "doubleThreshold", "fixed")) {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    signal <- as.vector(signal)
    Phases <- 2 * pi * (seq_len(phaseCount) - 1) / phaseCount
    Frequencies <- as.double(initialFrequencies)
    Schedule <- function(k, residue) {
        if (k > length(Frequencies)) return(NULL)
        list(cyclesPerSample = Frequencies[k] * samplingInterval,
             amplitude = .computeMaskAmplitude(residue, rule = "deviation",
                                               amplitudeRatio = amplitudeRatio))
    }
    Updates <- 0L
    repeat {
        Inner <- .decomposeMasked(signal, schedule = Schedule, phases = Phases, sifting = sifting,
                                   acceptIncomplete = TRUE)
        if (ncol(Inner$modes) < length(Frequencies)) {
            stop(errorCondition(paste("the masked decomposition has", ncol(Inner$modes),
                                      "modes for", length(Frequencies),
                                      "masks, so the mask update is undefined"),
                                class = c("xemdInsufficientModes", "moduleFailure")))
        }
        Attributes <- computeInstantaneousAttributes(Inner$modes,
                                                     samplingInterval = samplingInterval)
        Defined <- !is.na(Attributes$frequency)
        Power <- Attributes$amplitude^2 * Defined
        if (any(colSums(Power) == 0)) {
            stop(errorCondition(paste("a mode has zero squared amplitude where its frequency is",
                                      "defined, so its mask update is undefined"),
                                class = c("xemdZeroEnergy", "moduleFailure")))
        }
        Updated <- colSums(Power * replace(Attributes$frequency, !Defined, 0)) / colSums(Power)
        Updates <- Updates + 1L
        if (max(abs(Updated - Frequencies) / Frequencies) < tolerance) {
            StopReason <- "tolerance"
            break
        }
        if (Updates == maxIterations) {
            StopReason <- "iterationLimit"
            break
        }
        if (any(Updated <= 0)) {
            stop("the updated mask frequency of mode ", which(Updated <= 0)[1L],
                 " is not positive, so its relative change is undefined")
        }
        Frequencies <- Updated
    }
    list(modes = Inner$modes, residue = Inner$residue, iterations = Inner$iterations,
         converged = Inner$converged, stopReason = StopReason, maskFrequencies = Frequencies,
         maskIterations = Updates, innerStopReason = Inner$stopReason)
}
