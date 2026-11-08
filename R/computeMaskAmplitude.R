# Amplitude of the masks of one stage, in the units of the residue: under
# "maxima", amplitudeRatio V_m times the mean absolute value of the maxima of
# the residue (eq-mascaras-diadas-amplitud, mascaras.model:83-89); under
# "deviation", amplitudeRatio epsilon_0 times the standard deviation of the
# residue with divisor N (eq-mascaras-amplitud, mascaras.model:158-164).
.computeMaskAmplitude <- function(residue, rule, amplitudeRatio) {
    if (rule == "maxima") {
        return(amplitudeRatio * mean(abs(.findExtrema(residue)$maxima[, "value"])))
    }
    if (rule == "deviation") return(amplitudeRatio * sqrt(mean((residue - mean(residue))^2)))
    stop("unrecognized mask amplitude rule: ", rule)
}
