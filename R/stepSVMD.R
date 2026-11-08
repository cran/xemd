# One iteration of the extraction of a mode in SVMD
# (vmd.extensiones.es.qmd:257-314) on the full grid of |omega| in rad/s, with
# the spectrum of x, the mode spectrum, its center in rad/s and the multiplier
# of the current state, the current alpha in s^2, the dual step tau and the
# penalty of the previous centers. With A = alpha^2 (omega - center)^4 at the
# current center, the mode follows eq-svmd-modo-eliminado, whose denominator is
# the product of 1 + A and 1 + 2 alpha (omega - center)^2 + penalty; its
# centroid (eq-svmd-centro-aproximado) is the new center; the rest
# (eq-svmd-resto-actualizado) takes A at that new center; and the dual ascent
# (eq-svmd-dual) closes the iteration.
.stepSVMD <- function(spectrum, modeSpectrum, center, multiplier,
                      angularFrequency, alpha, tau, penalty) {
    Distance <- (angularFrequency - center)^2
    Factor <- (alpha * Distance)^2
    modeSpectrum <- (spectrum + Factor * modeSpectrum + multiplier / 2) /
        ((1 + Factor) * (1 + 2 * alpha * Distance + penalty))
    center <- .computeCenter(power = Mod(modeSpectrum)^2,
                             angularFrequency = angularFrequency)
    Factor <- (alpha * (angularFrequency - center)^2)^2
    Rest <- Factor / (1 + Factor) * (spectrum - modeSpectrum + multiplier / 2)
    list(modeSpectrum = modeSpectrum, center = center,
         multiplier = multiplier + tau * (spectrum - modeSpectrum - Rest))
}
