# Modes before the largest drop of the decreasing energies
# (eq-seleccion-potencia, eq-seleccion-orden, eq-seleccion-salto), most
# energetic first. The source fixes no tie-break, so a largest drop that is
# zero or not unique stops.
selectModesEnergy <- function(modes) {
    if (!is.matrix(modes) || !is.numeric(modes) || ncol(modes) < 2L ||
        !all(is.finite(modes))) {
        stop("'modes' must be a finite numeric matrix with at least two columns")
    }
    attributes(modes) <- list(dim = dim(modes))
    Energy <- colSums(modes^2)
    Order <- order(Energy, decreasing = TRUE)
    Drop <- -diff(Energy[Order])
    if (max(Drop) == 0) stop("every energy drop is zero")
    if (sum(Drop == max(Drop)) > 1L) stop("the largest energy drop is not unique")
    Order[seq_len(which.max(Drop))]
}
