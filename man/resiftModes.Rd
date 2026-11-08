\name{resiftModes}
\alias{resiftModes}
\encoding{UTF-8}
\title{Consecutive Post-Processing of Ensemble Modes}
\description{
Re-sifts a sequence of modes, typically the averaged modes of an ensemble,
by extracting the first sifted component of consecutive sums and
carrying each remainder to the next mode.
}
\usage{
resiftModes(modes, sifting)
}
\arguments{
\item{modes}{Finite numeric matrix, samples by modes, with at least two
columns \eqn{D_1, \ldots, D_K}{D_1, ..., D_K}.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected. The
chain has \eqn{K - 1} extractions, so \code{sifting$maxModes} does not
apply.}
}
\details{
With \eqn{u_1 = D_1 + D_2}, stage \eqn{i} extracts
\eqn{a_i = E_1(u_i)}, keeps the remainder \eqn{q_i = u_i - a_i} and adds
the next mode, \eqn{u_{i+1} = q_i + D_{i+2}}. The completed chain returns
\eqn{K - 1} modes and the residue \eqn{q_{K-1}} with
\code{"chainComplete"}, and
\deqn{\sum_{i=1}^{K-1} a_i + q_{K-1} = \sum_{j=1}^{K} D_j.}{sum_i a_i + q_{K-1} = sum_j D_j.}
When \eqn{E_1(u_i)} does not exist, the chain ends with
\code{"insufficientExtrema"}, the modes \eqn{a_1, \ldots, a_{i-1}}{a_1, ..., a_{i-1}} and the
residue \eqn{u_i + \sum_{j > i + 1} D_j}, which keeps the same sum.
}
\value{
A list with
\item{modes}{Double matrix, samples by the extracted modes.}
\item{residue}{Numeric vector, the remainder of the chain.}
\item{iterations}{Integer vector, the subtractions of each extraction.}
\item{converged}{Logical vector by extraction, as in
\code{\link{decomposeEMD}}.}
\item{stopReason}{\code{"chainComplete"} or \code{"insufficientExtrema"}.}
}
\references{
The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeEEMD}}, \code{\link{decomposeEMD}}
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
Modes <- decomposeEEMD(Signal, noise = Noise[, 1:2], amplitudeRatio = 0.2,
                       divisor = "realizations", sifting = Sifting)$modes
OUT <- resiftModes(Modes, sifting = Sifting)
c(input = ncol(Modes), output = ncol(OUT$modes))
OUT$stopReason
max(abs(rowSums(OUT$modes) + OUT$residue - rowSums(Modes)))
}
