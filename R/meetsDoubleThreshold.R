# Double threshold on sigma_h (emd.model:183-212): the fraction of the N
# samples with sigma_h >= theta_1 is at most alpha and no sample has
# sigma_h >= theta_2. Undefined samples (NA) are no exceedance and stay in the
# denominator N. The count condition |N_e - N_z| <= 1 belongs to the caller.
.meetsDoubleThreshold <- function(deviation, thresholds, fraction) {
    sum(deviation >= thresholds[1L], na.rm = TRUE) / length(deviation) <= fraction &&
        !any(deviation >= thresholds[2L], na.rm = TRUE)
}
