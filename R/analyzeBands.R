# Analysis coefficients of EWT, eq-ewt-analisis (ewt.sst.model.es.qmd:97-101):
# W_j = IDFT(conj(H_j) X) on the extension for each column of the bank
# (descomposicion.model.es.qmd:123-127). The bank is real and even in omega
# and the extension real, so W_j is real up to rounding.
.analyzeBands <- function(spectrum, filters) {
    Re(stats::mvfft(filters * spectrum, inverse = TRUE)) / nrow(filters)
}
