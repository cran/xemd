# Power-weighted centroid of eq-vmd-centro and eq-mvmd-centro on the full
# angular grid, in rad/s. Summing over the T bins with |omega| equals the
# half-line integrals with half weight at the zero and Nyquist bins. A mode
# without energy leaves the centroid undefined (vmd.model.es.qmd:209).
.computeCenter <- function(power, angularFrequency) {
    Energy <- sum(power)
    if (Energy == 0) {
        stop(errorCondition(
            "zero spectral energy: the center frequency is undefined",
            class = c("xemdZeroEnergy", "moduleFailure")))
    }
    sum(angularFrequency * power) / Energy
}
