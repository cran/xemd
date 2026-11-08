\name{decomposeIteratedMaskedEMD}
\alias{decomposeIteratedMaskedEMD}
\encoding{UTF-8}
\title{Iterated Masking Empirical Mode Decomposition}
\description{
Iterates a uniform-phase masked decomposition, updating each mask frequency
to the amplitude-weighted mean instantaneous frequency of its mode until a
relative-change test passes or the iteration limit is reached (itEMD).
}
\usage{
decomposeIteratedMaskedEMD(signal, samplingInterval, initialFrequencies,
                           phaseCount, amplitudeRatio, tolerance = 0.1,
                           maxIterations = 15, sifting)
}
\arguments{
\item{signal}{Numeric vector of finite real values without \code{dim}, in
the units of the quantity.}
\item{samplingInterval}{Sampling interval in seconds, a finite number
greater than zero, without default.}
\item{initialFrequencies}{Numeric vector of \eqn{K \ge 1} initial mask
frequencies in Hz, each in \eqn{(0, 1/(2\Delta t)]}, without default; its
length is the number of masks.}
\item{phaseCount}{Whole number \eqn{P \ge 2} of mask phases, without
default; the package tests use 4.}
\item{amplitudeRatio}{Dimensionless \eqn{\varepsilon_0 > 0}{epsilon_0 > 0}, without
default; the mask amplitude of each stage is \eqn{\varepsilon_0}{epsilon_0} times the
standard deviation of the residue with divisor \eqn{N}. The package
tests use 1.}
\item{tolerance}{Threshold of the largest relative change of the mask
frequencies, a finite number greater than zero; default \code{0.1}.}
\item{maxIterations}{Whole number \eqn{\ge 1} of frequency updates;
default \code{15}. The defaults of \code{tolerance} and
\code{maxIterations} are those of \code{iterated_mask_sift} in emd 0.8.1,
the code of a coauthor of itEMD.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected.}
}
\details{
Iteration \eqn{q} decomposes the signal with the uniform-phase masks of
\code{\link{decomposeMaskedEMD}}, one stage per frequency of
\eqn{\mathbf{f}^{(q)}}{f^{(q)}}, phases \eqn{2 \pi p / P} and the amplitude above.
Each frequency is then updated to the mean instantaneous frequency of its
mode weighted by the squared instantaneous amplitude, over the samples where
the frequency is defined,
\deqn{f_j^{(q+1)} = \frac{\sum_n A_j[n]^2 \, \nu_j[n]}{\sum_n A_j[n]^2},}{f_j^(q+1) = sum(A_j^2 nu_j) / sum(A_j^2),}
with \eqn{A_j} and \eqn{\nu_j} from
\code{\link{computeInstantaneousAttributes}}. The iteration stops with
\code{"tolerance"} when every relative change
\eqn{|f_j^{(q+1)} - f_j^{(q)}| / f_j^{(q)}} is below \code{tolerance}, or with
\code{"iterationLimit"} after \code{maxIterations} updates. The result is the
last decomposition together with the frequencies that produced it.
}
\value{
A list with the fields of \code{\link{decomposeMaskedEMD}} for the last
decomposition (\code{modes}, \code{residue}, \code{iterations},
\code{converged}) and
\item{stopReason}{\code{"tolerance"} or \code{"iterationLimit"}.}
\item{maskFrequencies}{Numeric vector, the \eqn{K} mask frequencies in Hz
that produced the returned decomposition.}
\item{maskIterations}{Integer, the number of frequency updates.}
\item{innerStopReason}{The termination of the returned masked
decomposition, with the vocabulary of \code{\link{decomposeMaskedEMD}}.}
}
\section{Conditions}{
When the masked decomposition has fewer modes than masks the update is
undefined and the function stops with an error of class
\code{xemdInsufficientModes}. A mode with zero squared amplitude where its
frequency is defined stops with an error of class \code{xemdZeroEnergy}. An
updated frequency that is not positive stops with an error when the
iteration continues.
}
\references{
Fabus, M. S., Quinn, A. J., Warnaby, C. E. and Woolrich, M. W. (2021).
Automatic decomposition of electrophysiological data into distinct
nonsinusoidal oscillatory modes. \emph{Journal of Neurophysiology},
126(5), 1670--1684. \doi{10.1152/jn.00315.2021}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeMaskedEMD}}, \code{\link{computeInstantaneousAttributes}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = 3)
# Three dyadic initial frequencies, with one sample per second.
OUT <- decomposeIteratedMaskedEMD(Signal, samplingInterval = 1,
                                  initialFrequencies = 55 / 3000 / 2 / 2^(0:2),
                                  phaseCount = 4, amplitudeRatio = 1,
                                  tolerance = 0.1, maxIterations = 15,
                                  sifting = Sifting)
OUT$maskFrequencies
OUT$maskIterations
OUT$stopReason
}
