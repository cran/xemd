\name{decomposeCEEMDAN}
\alias{decomposeCEEMDAN}
\encoding{UTF-8}
\title{Complete Ensemble Empirical Mode Decomposition with Adaptive Noise}
\description{
Extracts each mode as the average first mode of the common residue perturbed
by the modes of supplied noise realizations, with a constant ratio of the
noise to the residue at each stage.
}
\usage{
decomposeCEEMDAN(signal, noise, amplitudeRatio = 0.2, sifting)
}
\arguments{
\item{signal}{Numeric vector of at least two finite values without
\code{dim}, in the units of the quantity.}
\item{noise}{Finite numeric matrix with one row per sample of \code{signal}
and one column per realization, \eqn{M \ge 1}; column \eqn{m} is a
realization \eqn{w^{(m)}} of unit variance, used as supplied. The function
draws no random numbers.}
\item{amplitudeRatio}{Dimensionless ratio \eqn{\varepsilon \ge 0}{epsilon >= 0} of the
perturbation to the standard deviation of the residue at each stage
(divisor \eqn{N - 1}). Default \code{0.2}, as in
\code{\link{decomposeEEMD}}.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected. The
same list decomposes the noise realizations.}
}
\details{
Each column of \code{noise} is decomposed once with \code{sifting},
\eqn{\eta_k^{(m)} = E_k(w^{(m)})}, and \eqn{\eta_0^{(m)} = w^{(m)}}. From
\eqn{r_0 = x}, before stage \eqn{k} the common residue \eqn{r_{k-1}} ends
the decomposition, in this order, with \code{"terminalResidue"} when it has
at most one extremum, with \code{"modeLimit"} when \eqn{k - 1} equals
\code{sifting$maxModes} and with \code{"insufficientExtrema"} when its
envelopes cannot be built. Stage \eqn{k} perturbs \eqn{r_{k-1}} with
\eqn{p^{(m)} = \beta_{k-1}^{(m)} \eta_{k-1}^{(m)}},
\eqn{\beta_{k-1}^{(m)} = \varepsilon \, \mathrm{sd}(r_{k-1}) / \mathrm{sd}(\eta_{k-1}^{(m)})}{beta_{k-1}^{(m)} = epsilon sd(r_{k-1}) / sd(eta_{k-1}^{(m)})},
and computes
\deqn{d_k = \frac{1}{M} \sum_{m=1}^{M} E_1\left(r_{k-1} + p^{(m)}\right), \qquad r_k = r_{k-1} - d_k.}{d_k = (1/M) sum_m E_1(r_{k-1} + p^(m)), r_k = r_{k-1} - d_k.}
A realization without \eqn{\eta_{k-1}} contributes a null perturbation, and
a perturbed residue without a first mode contributes zero; both stay in the
divisor \eqn{M}.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes.}
\item{residue}{Numeric vector, the final common residue.}
\item{iterations}{Integer vector by mode, the sum of the subtractions of
the \eqn{M} extractions of the stage.}
\item{converged}{Logical vector by mode: \code{TRUE} when the \eqn{M}
extractions accepted, \code{FALSE} when one ended without accepting or
could not begin, \code{NA} when all completed fixed sifting.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}.}
}
\section{Conditions}{
A noise mode of zero standard deviation leaves its scale undefined and
stops the function with an error of class \code{xemdZeroEnergy}.
}
\references{
Torres, M. E., Colominas, M. A., Schlotthauer, G. and Flandrin, P. (2011).
A complete ensemble empirical mode decomposition with adaptive noise. In
\emph{2011 IEEE International Conference on Acoustics, Speech and Signal
Processing (ICASSP)}, 4144--4147. \doi{10.1109/ICASSP.2011.5947265}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeICEEMDAN}}, \code{\link{decomposeEEMD}},
\code{\link{decomposeCEEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 1000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
set.seed(628)
Noise <- t(matrix(stats::rnorm(10L * length(Signal)), nrow = 10L,
                  byrow = TRUE))
# Two of the ten realizations and three modes keep the example short.
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = 3)
OUT <- decomposeCEEMDAN(Signal, noise = Noise[, 1:2], amplitudeRatio = 0.2,
                        sifting = Sifting)
dim(OUT$modes)
OUT$stopReason
max(abs(rowSums(OUT$modes) + OUT$residue - Signal))
}
