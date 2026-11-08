# Inverse of .transformReflection for each column: inverse DFT divided by T
# (eq-marco-idft) at the observed positions of the extension
# (vmd.model.es.qmd:269). The spectra are Hermitian, so the discarded
# imaginary part is rounding.
.restoreObserved <- function(spectrum, observed) {
    Re(stats::mvfft(spectrum, inverse = TRUE)[observed, , drop = FALSE]) /
        nrow(spectrum)
}
