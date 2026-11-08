# Empirical mode decomposition by successive extraction (emd.model:289-337;
# descomposicion.model:143-157). r_0 is the signal; before extraction k the
# decomposition ends with "terminalResidue" when r_{k-1} has at most one
# extremum, with "modeLimit" when k - 1 equals maxModes, and with
# "insufficientExtrema" when E_1 is undefined on r_{k-1}; otherwise
# c_k = E_1(r_{k-1}) and r_k = r_{k-1} - c_k. The signal is a vector of
# samples in the units of the quantity; positions are sample indices. The
# decomposition runs on the bare samples, so no attribute of signal (names,
# tsp, class) reaches the outputs. sifting is completed and validated by
# .checkSifting; "consecutive" belongs to the vector sifting and is rejected.
decomposeEMD <- function(signal, sifting) {
    if (!is.numeric(signal) || !is.null(dim(signal)) || !all(is.finite(signal))) {
        stop("signal must be a numeric vector of finite real values without dim")
    }
    signal <- as.vector(signal)
    sifting <- .checkSifting(sifting)
    if (sifting$stopRule == "consecutive") {
        stop("sifting$stopRule must be SD, Snumber, doubleThreshold or fixed for scalar sifting")
    }
    .decomposeScalar(signal, sifting = sifting, acceptIncomplete = TRUE)
}
