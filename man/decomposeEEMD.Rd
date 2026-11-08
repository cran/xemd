\name{decomposeEEMD}
\alias{decomposeEEMD}
\encoding{UTF-8}
\title{Ensemble Empirical Mode Decomposition}
\description{
Decomposes noise-perturbed copies of a signal by scalar sifting and averages
their completed modes by extraction order.
}
\usage{
decomposeEEMD(signal, noise, amplitudeRatio = 0.2, divisor, sifting)
}
\arguments{
\item{signal}{Numeric vector of at least two finite values without
\code{dim}, in the units of the quantity.}
\item{noise}{Finite numeric matrix with one row per sample of \code{signal}
and one column per realization, \eqn{M \ge 1}. Column \eqn{m} is a
realization \eqn{z^{(m)}} of unit variance, used as supplied: the function
neither generates nor normalizes it, and it draws no random numbers.}
\item{amplitudeRatio}{Dimensionless ratio \eqn{\varepsilon \ge 0}{epsilon >= 0} of the
noise standard deviation to the standard deviation of \code{signal}
(divisor \eqn{N - 1}). Default \code{0.2}, the value common to libeemd and
emd 0.8.1.}
\item{divisor}{Categorical, without default: \code{"realizations"} divides
the sum of the modes of order \eqn{j} by \eqn{M}; \code{"available"}
divides it by the number \eqn{M_j} of realizations that have that mode.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected.}
}
\details{
Realization \eqn{m} is \eqn{x + w^{(m)}} with
\eqn{w^{(m)} = \varepsilon \, \mathrm{sd}(x) \, z^{(m)}}{w^{(m)} = epsilon sd(x) z^{(m)}}; it is decomposed
by scalar sifting with \code{sifting}, in the order of the columns of
\code{noise}. Each realization stops before adding an extraction whose sifting lost
its envelopes and keeps the input of that attempt in its residue. Accepted
adaptive extractions, completed fixed sifting and extractions that reached
the iteration limit are retained, the last ones flagged \code{FALSE}. Absent modes are extended with zeros, and
mode \eqn{j} is
\deqn{D_j = \frac{1}{M} \sum_{m=1}^{M} c_j^{(m)}}{D_j = (1/M) sum_m c_j^(m)}
with \code{"realizations"}, or the same sum divided by \eqn{M_j} with
\code{"available"}. The residue is the mean of the \eqn{M} final residues
and \code{meanNoise} the mean of the added perturbations, so that with
\code{"realizations"} the modes, the residue and \code{signal} satisfy
\code{rowSums(modes) + residue = signal + meanNoise} up to rounding. A
realization without modes contributes zeros to the modes and its input to
the residue. Matching by extraction order does not establish a physical
correspondence between the modes of different realizations.

The method has no outer termination, so the result has no
\code{stopReason}. The number of realizations has no default; the audited
libraries use 100 (PyEMD), 250 (libeemd) and 4 (emd 0.8.1). hht divides by
\eqn{M} and PyEMD by \eqn{M_j}.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes, with as many columns as the
largest number of modes among the realizations.}
\item{residue}{Numeric vector, the mean final residue.}
\item{meanNoise}{Numeric vector, the mean added perturbation, in the units
of \code{signal}.}
\item{iterations}{Integer vector by mode, the sum of the subtractions of
the realizations that have the mode.}
\item{converged}{Logical vector by mode: \code{TRUE} only when every
realization that has the mode accepted it; \code{NA} for completed fixed
sifting. Failed attempts are reported in \code{termination}.}
\item{termination}{Data frame with one row per realization:
\code{realization}, its column index in \code{noise}; \code{mode}, the next
extraction index after the retained modes; \code{iterations}, subtractions
of the interrupted attempt, zero if it did not start; and \code{stopReason},
\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}. These rows describe individual decompositions,
not the averaged mode indices.}
}
\references{
Wu, Z. and Huang, N. E. (2009). Ensemble empirical mode decomposition: A
noise-assisted data analysis method. \emph{Advances in Adaptive Data
Analysis}, 1(1), 1--41. \doi{10.1142/S1793536909000047}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeCEEMD}}, \code{\link{decomposeCEEMDAN}},
\code{\link{decomposeICEEMDAN}}, \code{\link{resiftModes}},
\code{\link{measureDispersion}}
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
OUT <- decomposeEEMD(Signal, noise = Noise[, 1:2], amplitudeRatio = 0.2,
                     divisor = "realizations", sifting = Sifting)
dim(OUT$modes)
OUT$iterations
max(abs(rowSums(OUT$modes) + OUT$residue - Signal - OUT$meanNoise))
}
