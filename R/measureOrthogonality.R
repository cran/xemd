# Inner products between the columns (eq-validacion-energia-cruzada), without
# the sampling interval of eq-marco-producto-discreto, which cancels in the
# cosine. The cosine keeps its sign.
measureOrthogonality <- function(components) {
    if (!is.matrix(components) || !is.numeric(components) ||
        !all(is.finite(components))) {
        stop("'components' must be a finite numeric matrix")
    }
    attributes(components) <- list(dim = dim(components))
    InnerProduct <- crossprod(components)
    Norm <- sqrt(diag(InnerProduct))
    if (any(Norm == 0)) {
        stop(errorCondition("'components' has a column of zero norm",
                            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    list(innerProduct = InnerProduct, cosine = InnerProduct / outer(Norm, Norm))
}
