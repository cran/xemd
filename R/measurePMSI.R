# Pseudo-mode splitting index of each pair of adjacent modes
# (eq-validacion-pmsi). The regularized denominator adds the 1e-8 of the UPEMD
# definition, an absolute scale absent from the homogeneous form.
measurePMSI <- function(modes, denominator) {
    if (!is.matrix(modes) || !is.numeric(modes) || ncol(modes) < 1L ||
        !all(is.finite(modes))) {
        stop("'modes' must be a finite numeric matrix with at least one column")
    }
    if (!identical(denominator, "homogeneous") && !identical(denominator, "regularized")) {
        stop("'denominator' must be \"homogeneous\" or \"regularized\"")
    }
    attributes(modes) <- list(dim = dim(modes))
    k <- ncol(modes)
    Energy <- colSums(modes^2)
    Denominator <- Energy[-k] + Energy[-1L]
    if (denominator == "regularized") Denominator <- Denominator + 1e-8
    if (any(Denominator == 0)) {
        stop(errorCondition("a pair of adjacent modes has zero energy",
                            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    pmax(colSums(modes[, -k, drop = FALSE] * modes[, -1L, drop = FALSE]) / Denominator, 0)
}
