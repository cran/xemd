\name{decomposeCEEMD}
\alias{decomposeCEEMD}
\encoding{UTF-8}
\title{Complementary Ensemble Empirical Mode Decomposition}
\description{
Decomposes pairs of copies of a signal perturbed with the same noise
realization of opposite signs and averages the modes of all the copies by
extraction order.
}
\usage{
decomposeCEEMD(signal, noise, amplitudeRatio = 0.2, sifting)
}
\arguments{
\item{signal}{Numeric vector of at least two finite values without
\code{dim}, in the units of the quantity.}
\item{noise}{Finite numeric matrix with one row per sample of \code{signal}
and one column per realization, \eqn{M \ge 1}; column \eqn{m} is a
realization \eqn{z^{(m)}} of unit variance, used as supplied. The function
draws no random numbers.}
\item{amplitudeRatio}{Dimensionless ratio \eqn{\varepsilon \ge 0}{epsilon >= 0}; the
perturbation is \eqn{a z^{(m)}} with \eqn{a = \varepsilon \, \mathrm{sd}(x)}{a = epsilon sd(x)}
(divisor \eqn{N - 1}). Default \code{0.2}, as in
\code{\link{decomposeEEMD}}.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected.}
}
\details{
For each column, the inputs \eqn{x + a z^{(m)}} and \eqn{x - a z^{(m)}} are
decomposed by scalar sifting, the positive input first. Each stops before
adding an extraction whose sifting lost its envelopes and keeps its input
in the residue. Accepted adaptive extractions, completed fixed sifting and
extractions that reached the iteration limit are retained.
The \eqn{2M} decompositions are extended with zero modes to the
largest number of modes and averaged by extraction order,
\deqn{d_j = \frac{1}{2M} \sum_{m=1}^{M} \left(c_j^{(m,+)} + c_j^{(m,-)}\right),}{d_j = (1/(2M)) sum_m (c_j^(m,+) + c_j^(m,-)),}
and the residue is the mean of the \eqn{2M} residues, so the perturbations
cancel in the reconstruction: \code{rowSums(modes) + residue} equals
\code{signal} up to rounding. The method has no outer termination and
returns neither \code{meanNoise} nor \code{stopReason}.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes.}
\item{residue}{Numeric vector, the mean residue of the \eqn{2M}
decompositions.}
\item{iterations}{Integer vector by mode, the sum of the subtractions of
the decompositions that have the mode.}
\item{converged}{Logical vector by mode: \code{TRUE} only when every
decomposition that has the mode accepted it; \code{NA} for completed fixed
sifting. Failed attempts are reported in \code{termination}.}
\item{termination}{Data frame with one row per signed realization:
\code{realization}, its column index in \code{noise}; \code{sign}, 1 or -1;
and \code{mode}, \code{iterations} and \code{stopReason} as in
\code{\link{decomposeEEMD}}. A stopped realization stays in its pair and
in the divisor \eqn{2M}; its remaining content stays in the mean residue.}
}
\references{
Yeh, J.-R., Shieh, J.-S. and Huang, N. E. (2010). Complementary ensemble
empirical mode decomposition: A novel noise enhanced data analysis method.
\emph{Advances in Adaptive Data Analysis}, 2(2), 135--156.
\doi{10.1142/S1793536910000422}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeEEMD}}, \code{\link{decomposeCEEMDAN}},
\code{\link{decomposeICEEMDAN}}
}
\examples{
Time <- seq(0, 9, length.out = 1000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
set.seed(628)
Noise <- t(matrix(stats::rnorm(10L * length(Signal)), nrow = 10L,
                  byrow = TRUE))
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = Inf)
# Two of the ten realizations keep the example short.
OUT <- decomposeCEEMD(Signal, noise = Noise[, 1:2], amplitudeRatio = 0.2,
                      sifting = Sifting)
dim(OUT$modes)
max(abs(rowSums(OUT$modes) + OUT$residue - Signal))
}
