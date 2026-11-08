# Deviation sigma_h = |m / a| of the envelope mean m relative to the envelope
# amplitude a at every sample (emd.model:169-181); the caller forms m and a.
# Where a = 0 the ratio is undefined when m = 0, returned as NA, and +Inf when
# m != 0, as in emd.m and emd 0.8.1 (emd.rilling.implementation:189;
# sift.implementation:5).
.measureDeviation <- function(mean, amplitude) {
    OUT <- abs(mean / amplitude)
    OUT[mean == 0 & amplitude == 0] <- NA_real_
    OUT
}
