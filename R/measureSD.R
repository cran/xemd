# Criterion SD of the 1998 formulation (emd.model:106-124): the sum over all N
# samples of (h_{j-1} - h_j)^2 / h_{j-1}^2, without the factor 1/N. A null
# denominator leaves an infinite or undefined term; it is not regularized.
.measureSD <- function(previous, current) {
    sum((previous - current)^2 / previous^2)
}
