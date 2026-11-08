# Empirical wavelet transform (ewt.sst.model.es.qmd:7-123) of a signal sampled
# every samplingInterval seconds into at most bandCount bands: boundaries at
# the midpoints between the dominant maxima of its Fourier magnitude, fewer
# bands when fewer maxima exist; the filter bank on the mirror extension; the
# analysis coefficients and, through the dual bank, the components c_j on the
# input grid, which add up to the signal and are not IMF. boundaries are in Hz.
# A detection list selects the boundary detection of the EWT code of Gilles
# (supplement ewt): detector "locmax", "locmaxmin" or "locmaxminf", required;
# logarithm, FALSE by default; regularization "none", "average" or "gaussian",
# "gaussian" by default, with window 10 and sigma 1.5 by default, the initial
# values of that code; and, with bandCount = NULL, thresholdRatio in (0, 1)
# for the automatic count, which needs locmax, no logarithm and no
# regularization. detection = NULL keeps the midpoint route.
decomposeEWT <- function(signal, samplingInterval, bandCount, detection = NULL) {
    stopifnot(
        "signal must be a numeric vector of at least four finite values" =
            is.numeric(signal) && is.null(dim(signal)) && length(signal) >= 4L &&
            all(is.finite(signal)),
        "samplingInterval must be a finite number greater than zero" =
            is.numeric(samplingInterval) && length(samplingInterval) == 1L &&
            is.finite(samplingInterval) && samplingInterval > 0,
        "bandCount must be NULL or a whole number not below two" =
            is.null(bandCount) ||
            (is.numeric(bandCount) && length(bandCount) == 1L &&
             is.finite(bandCount) && bandCount >= 2 &&
             bandCount == floor(bandCount)),
        "detection must be NULL or a named list" =
            is.null(detection) || (is.list(detection) && !is.null(names(detection))),
        "bandCount = NULL needs a detection list with thresholdRatio" =
            !is.null(bandCount) || !is.null(detection$thresholdRatio))
    if (!is.null(detection)) {
        FIELDS <- c("detector", "logarithm", "regularization", "window", "sigma",
                    "thresholdRatio")
        if (anyDuplicated(names(detection)) || !all(names(detection) %in% FIELDS)) {
            stop("detection must have distinct fields among ", paste(FIELDS, collapse = ", "))
        }
        AUX <- list(logarithm = FALSE, regularization = "gaussian", window = 10,
                    sigma = 1.5)
        AUX <- AUX[setdiff(names(AUX), names(detection))]
        detection[names(AUX)] <- AUX
        isNumber <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
        stopifnot(
            "detection$detector must be \"locmax\", \"locmaxmin\" or \"locmaxminf\"" =
                is.character(detection$detector) && length(detection$detector) == 1L &&
                detection$detector %in% c("locmax", "locmaxmin", "locmaxminf"),
            "detection$logarithm must be TRUE or FALSE" =
                isTRUE(detection$logarithm) || isFALSE(detection$logarithm),
            "detection$regularization must be \"none\", \"average\" or \"gaussian\"" =
                is.character(detection$regularization) &&
                length(detection$regularization) == 1L &&
                detection$regularization %in% c("none", "average", "gaussian"),
            "detection$window must be a whole number not below one" =
                isNumber(detection$window) && detection$window >= 1 &&
                detection$window == floor(detection$window),
            "detection$sigma must be a finite number greater than zero" =
                isNumber(detection$sigma) && detection$sigma > 0,
            "detection$thresholdRatio must be a number between zero and one, exclusive" =
                is.null(detection$thresholdRatio) ||
                (isNumber(detection$thresholdRatio) && detection$thresholdRatio > 0 &&
                 detection$thresholdRatio < 1),
            "detection$thresholdRatio applies only with bandCount = NULL" =
                is.null(bandCount) || is.null(detection$thresholdRatio),
            "bandCount = NULL needs locmax without logarithm or regularization" =
                !is.null(bandCount) ||
                (detection$detector == "locmax" && !detection$logarithm &&
                 detection$regularization == "none"))
    }
    # The outputs inherit no attribute of signal (names, tsp, class).
    signal <- as.vector(signal)
    Boundaries <- .segmentSpectrum(signal = signal, bandCount = bandCount,
                                   detection = detection)
    Transform <- .transformEWT(signal)
    Filters <- .buildFilterBank(boundaries = Boundaries,
                                size = length(Transform$spectrum))
    Coefficients <- .analyzeBands(spectrum = Transform$spectrum, filters = Filters)
    list(modes = .synthesizeBands(coefficients = Coefficients, filters = Filters,
                                  observed = Transform$observed),
         boundaries = Boundaries / (2 * pi * samplingInterval))
}
