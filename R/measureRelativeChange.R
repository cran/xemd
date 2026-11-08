# Relative change between consecutive iterates of iterative filtering,
# SD_q = ||f^(q) - f^(q-1)||_2 / ||f^(q-1)||_2 (eq-if-parada-relativa;
# lmd.model:290-297): a ratio of norms, not of squared norms. By Parseval it
# takes the same value on two real states or on their full discrete spectra.
# A null previous state leaves it undefined (NaN).
.measureRelativeChange <- function(previous, current) {
    sqrt(sum(Mod(current - previous)^2)) / sqrt(sum(Mod(previous)^2))
}
