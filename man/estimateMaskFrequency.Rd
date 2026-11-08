\name{estimateMaskFrequency}
\alias{estimateMaskFrequency}
\encoding{UTF-8}
\title{Mask Frequency of the 2005 Masking}
\description{
Estimates the frequency of the masking signal of 2005 from the
instantaneous amplitude and frequency of the first sifted component
of a signal.
}
\usage{
estimateMaskFrequency(signal, samplingInterval, sifting)
}
\arguments{
\item{signal}{Numeric vector of finite real values without \code{dim}, in
the units of the quantity.}
\item{samplingInterval}{Sampling interval in seconds, a finite number
greater than zero, without default.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected.}
}
\details{
With \eqn{A_1} and \eqn{\nu_1} the instantaneous amplitude and frequency,
computed by \code{\link{computeInstantaneousAttributes}}, of the first
sifted component \eqn{y_1 = E_1(x)} without mask, the estimate is
\deqn{\bar\nu = \frac{\sum_n A_1[n] \, \nu_1[n]^2}{\sum_n A_1[n] \, \nu_1[n]},}{nu = sum(A_1 nu_1^2) / sum(A_1 nu_1),}
with the sums over the samples where \eqn{\nu_1} is defined. The result is
the \code{frequency} of the rule \code{"pair"} of
\code{\link{decomposeMaskedEMD}}.
}
\value{
A number, the mask frequency in Hz.
}
\section{Conditions}{
A signal without a first extracted component stops with an error. When
the weights \eqn{A_1 \nu_1} sum to zero the estimate is undefined, and the
function stops with an error of class \code{xemdZeroEnergy}.
}
\references{
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
# One sample per second, so the result is in cycles per sample.
estimateMaskFrequency(Signal, samplingInterval = 1, sifting = Sifting)
# The time step of the signal, 9 / 2999 s.
estimateMaskFrequency(Signal, samplingInterval = 9 / 2999, sifting = Sifting)
}
