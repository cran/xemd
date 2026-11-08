# Masked empirical mode decomposition (mascaras.model:9-180) under one of
# three mask rules, each a named control of the method. "pair", the opposite
# sine masks of 2005: one extraction c_1 = (E_1(x + w) + E_1(x - w)) / 2 with
# w[n] = amplitude sin(2 pi f n dt), the phases -pi / 2 and pi / 2 of
# eq-mascaras-sinusoide (eq-mascaras-par, eq-mascaras-residuo,
# eq-mascaras-seno-original). "dyadic", EMD-DMS: phases 0 and pi, frequencies
# f_j = 2^(1 - j) f_1 and amplitude V_m times the mean absolute maximum of the
# residue (eq-mascaras-diadas, eq-mascaras-diadas-amplitud,
# eq-mascaras-diadas-extraccion). "uniformPhase", UPEMD: P phases 2 pi p / P,
# frequencies f_k = 1 / (2^k dt) and amplitude epsilon_0 sd_N of the residue
# (eq-mascaras-fases, eq-mascaras-periodo, eq-mascaras-amplitud,
# eq-mascaras-multinivel-modo). Frequencies are in Hz from samplingInterval in
# seconds. The schedule of "pair" has one stage, so it ends with "modeLimit"
# after one extraction unless the residue is terminal.
decomposeMaskedEMD <- function(signal, samplingInterval, masking, sifting) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal))) {
        stop("signal must be a numeric vector of finite real values without dim")
    }
    isPositive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
    if (!isPositive(samplingInterval)) {
        stop("samplingInterval must be a finite number of seconds greater than zero")
    }
    FIELDS <- list(pair = c("frequency", "amplitude"),
                   dyadic = c("firstFrequency", "amplitudeRatio"),
                   uniformPhase = c("amplitudeRatio", "phaseCount"))
    if (!is.list(masking)) stop("masking must be a named list")
    Rule <- masking[["rule"]]
    if (!(is.character(Rule) && length(Rule) == 1L && Rule %in% names(FIELDS))) {
        stop("masking$rule must be one of ", paste(names(FIELDS), collapse = ", "))
    }
    Required <- c("rule", FIELDS[[Rule]])
    if (is.null(names(masking)) || anyDuplicated(names(masking)) ||
        !setequal(names(masking), Required)) {
        stop("masking with rule ", Rule, " must have exactly the fields ",
             paste(Required, collapse = ", "))
    }
    isFrequency <- function(x) isPositive(x) && x <= 1 / (2 * samplingInterval)
    if (Rule == "pair" && !isFrequency(masking$frequency)) {
        stop("masking$frequency must be a frequency in Hz in (0, 1 / (2 samplingInterval)]")
    }
    if (Rule == "pair" && !isPositive(masking$amplitude)) {
        stop("masking$amplitude must be a finite number greater than zero")
    }
    if (Rule == "dyadic" && !isFrequency(masking$firstFrequency)) {
        stop("masking$firstFrequency must be a frequency in Hz in (0, 1 / (2 samplingInterval)]")
    }
    if (Rule %in% c("dyadic", "uniformPhase") && !isPositive(masking$amplitudeRatio)) {
        stop("masking$amplitudeRatio must be a finite number greater than zero")
    }
    if (Rule == "uniformPhase" &&
        !(isPositive(masking$phaseCount) && masking$phaseCount >= 2 &&
          masking$phaseCount == round(masking$phaseCount))) {
        stop("masking$phaseCount must be a whole number of at least two")
    }
    sifting <- .checkSifting(sifting)
    if (!sifting$stopRule %in% c("SD", "Snumber", "doubleThreshold", "fixed")) {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    if (Rule == "pair") {
        Phases <- c(-pi / 2, pi / 2)
        Schedule <- function(k, residue) {
            if (k > 1L) return(NULL)
            list(cyclesPerSample = masking$frequency * samplingInterval,
                 amplitude = masking$amplitude)
        }
    }
    if (Rule == "dyadic") {
        Phases <- c(0, pi)
        Schedule <- function(k, residue) {
            list(cyclesPerSample = 2^(1 - k) * masking$firstFrequency * samplingInterval,
                 amplitude = .computeMaskAmplitude(residue, rule = "maxima",
                                                   amplitudeRatio = masking$amplitudeRatio))
        }
    }
    if (Rule == "uniformPhase") {
        Phases <- 2 * pi * (seq_len(masking$phaseCount) - 1) / masking$phaseCount
        Schedule <- function(k, residue) {
            list(cyclesPerSample = 1 / 2^k,
                 amplitude = .computeMaskAmplitude(residue, rule = "deviation",
                                                   amplitudeRatio = masking$amplitudeRatio))
        }
    }
    OUT <- .decomposeMasked(as.vector(signal), schedule = Schedule, phases = Phases,
                            sifting = sifting, acceptIncomplete = FALSE)
    list(modes = OUT$modes, residue = OUT$residue, iterations = OUT$iterations,
         converged = OUT$converged, stopReason = OUT$stopReason,
         maskFrequencies = OUT$cyclesPerSample / samplingInterval, termination = OUT$termination)
}
