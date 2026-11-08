# Domain of the sifting configuration theta (eemd.model:43). The fields that
# have a default under the received stopRule are completed first; the
# configuration must then have exactly the fields of its rule, each in its
# domain (emd.model:106-226, :230-240; eemd.model:283-295). Counts are
# accepted as doubles with a whole value; under "Snumber" and "consecutive"
# the limit must allow the run of states that the rule requires, M >= S
# (emd.model:218). Returns the complete configuration; any departure is an
# error. Defaults: maxModes Inf extracts until the residue ends the
# decomposition; under "doubleThreshold" those of emd.m, reference code of the
# authors of the double threshold (emd.rilling.implementation:5, :57); under
# "fixed" reflection, common to PyEMD and emd 0.8.1 (pyemd.implementation:3,
# :187; sift.implementation:3, :31); under "consecutive" those of memd.m of
# Rehman and Mandic, the authors' code: n_iter = 2 (memd.m:40-44, :575-576),
# MAXITERATIONS = 1000 (:609-610) and reflection with two extrema per side
# (:187, :229).
.checkSifting <- function(sifting) {
    FIELDS <- list(SD = c("tolerance", "maxSift"), Snumber = c("stableSift", "maxSift"),
                   doubleThreshold = c("thresholds", "fraction", "maxSift"), fixed = "fixedSift",
                   consecutive = c("consecutiveCount", "maxSift"))
    DEFAULTS <- list(doubleThreshold = list(thresholds = c(0.05, 0.5), fraction = 0.05,
                                            maxSift = 2000, boundary = "reflection"),
                     fixed = list(boundary = "reflection"),
                     consecutive = list(consecutiveCount = 2, maxSift = 1000,
                                        boundary = "reflection"))
    if (!is.list(sifting)) stop("sifting must be a named list")
    Rule <- sifting[["stopRule"]]
    if (!(is.character(Rule) && length(Rule) == 1L && Rule %in% names(FIELDS))) {
        stop("sifting$stopRule must be one of ", paste(names(FIELDS), collapse = ", "))
    }
    AUX <- c(list(maxModes = Inf), DEFAULTS[[Rule]])
    AUX <- AUX[setdiff(names(AUX), names(sifting))]
    sifting[names(AUX)] <- AUX
    Required <- c("stopRule", FIELDS[[Rule]], "maxModes", "boundary")
    if (is.null(names(sifting)) || anyDuplicated(names(sifting)) || !setequal(names(sifting), Required)) {
        stop("sifting with stopRule ", Rule, " must have exactly the fields ",
             paste(Required, collapse = ", "))
    }
    isNumber <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    isCount <- function(x) isNumber(x) && x >= 1 && x == round(x)
    if (Rule == "SD" && !(isNumber(sifting$tolerance) && sifting$tolerance > 0)) {
        stop("sifting$tolerance must be a finite number greater than zero")
    }
    if (Rule == "Snumber" && !isCount(sifting$stableSift)) {
        stop("sifting$stableSift must be a whole number of at least one")
    }
    if (Rule == "consecutive" && !isCount(sifting$consecutiveCount)) {
        stop("sifting$consecutiveCount must be a whole number of at least one")
    }
    if (Rule == "doubleThreshold" &&
        !(is.numeric(sifting$thresholds) && length(sifting$thresholds) == 2L &&
          all(is.finite(sifting$thresholds)) && 0 < sifting$thresholds[1L] &&
          sifting$thresholds[1L] < sifting$thresholds[2L])) {
        stop("sifting$thresholds must be two finite numbers with 0 < thresholds[1] < thresholds[2]")
    }
    if (Rule == "doubleThreshold" &&
        !(isNumber(sifting$fraction) && sifting$fraction >= 0 && sifting$fraction < 1)) {
        stop("sifting$fraction must be a finite number in [0, 1)")
    }
    if (Rule == "fixed" && !isCount(sifting$fixedSift)) {
        stop("sifting$fixedSift must be a whole number of at least one")
    }
    if (Rule != "fixed" && !isCount(sifting$maxSift)) {
        stop("sifting$maxSift must be a whole number of at least one")
    }
    if (Rule == "Snumber" && sifting$maxSift < sifting$stableSift) {
        stop("sifting$maxSift must be at least sifting$stableSift")
    }
    if (Rule == "consecutive" && sifting$maxSift < sifting$consecutiveCount) {
        stop("sifting$maxSift must be at least sifting$consecutiveCount")
    }
    if (!(is.numeric(sifting$maxModes) && length(sifting$maxModes) == 1L &&
          (isTRUE(sifting$maxModes == Inf) || isCount(sifting$maxModes)))) {
        stop("sifting$maxModes must be a whole number of at least one or Inf")
    }
    if (!(is.character(sifting$boundary) && length(sifting$boundary) == 1L &&
          sifting$boundary %in% c("reflection", "extrapolation"))) {
        stop("sifting$boundary must be reflection or extrapolation")
    }
    sifting
}
