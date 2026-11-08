# Relative change between consecutive sweeps (eq-vmd-parada, eq-mvmd-parada):
# sum over modes and channels of ||U(n+1) - U(n)||^2 / ||U(n)||^2, with norms
# over the T bins of the T x K x C mode spectra. A coordinate without energy
# in the previous sweep leaves the ratio undefined; the change is then Inf, so
# the stopping test fails until the iteration limit.
.measureSweepChange <- function(current, previous) {
    Energy <- colSums(Mod(previous)^2)
    if (any(Energy == 0)) return(Inf)
    sum(colSums(Mod(current - previous)^2) / Energy)
}
