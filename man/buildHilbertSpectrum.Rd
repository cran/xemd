\name{buildHilbertSpectrum}
\alias{buildHilbertSpectrum}
\encoding{UTF-8}
\title{Hilbert Spectrum by Frequency Cells}
\description{
Accumulates the instantaneous amplitude, or its square, of every mode in the
frequency cell of its instantaneous frequency, sample by sample, and
integrates the result over time.
}
\usage{
buildHilbertSpectrum(frequency, amplitude, edges, exponent, samplingInterval)
}
\arguments{
\item{frequency}{Numeric matrix, samples by modes, of instantaneous
frequencies in Hz, finite or \code{NA}.}
\item{amplitude}{Finite non-negative matrix with the dimensions of
\code{frequency}, the instantaneous amplitudes in the units of the
signal.}
\item{edges}{At least two finite, strictly increasing cell edges in Hz;
cell \eqn{b} is \eqn{[e_b, e_{b+1})}. Negative values are admitted. No
default.}
\item{exponent}{The power \eqn{p}: 1 for the amplitude spectrum, 2 for the
spectrum weighted by squared amplitude. No default.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
}
\details{
The spectrum holds the mass of each cell,
\deqn{H[n, b] = \sum_{k:\, f_k[n] \in [e_b, e_{b+1})} a_k[n]^p,}{H[n, b] = sum over modes k with f_k[n] in [e_b, e_{b+1}) of a_k[n]^p,}
in units of the signal to the power \eqn{p}. A frequency outside the edges
or \code{NA} adds nothing. The marginal spectrum of a cell is the time
integral \eqn{\Delta t \sum_n H[n, b]}. No scaling or division by the cell
width is applied: the values are masses, not densities. The inputs are
usually the \code{frequency} and \code{amplitude} of
\code{\link{computeInstantaneousAttributes}}.
}
\value{
A list with
\item{spectrum}{Double matrix, samples by cells, in units of the signal to
the power \eqn{p}.}
\item{marginal}{Numeric vector by cell, in units of the signal to the power
\eqn{p} times seconds.}
}
\seealso{
\code{\link{computeInstantaneousAttributes}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Modes <- decomposeEMD(Signal,
                      sifting = list(stopRule = "doubleThreshold"))$modes
Attributes <- computeInstantaneousAttributes(Modes,
                                             samplingInterval = 9 / 2999)
OUT <- buildHilbertSpectrum(Attributes$frequency, Attributes$amplitude,
                            edges = c(1, 2, 3, 4, 5), exponent = 2,
                            samplingInterval = 9 / 2999)
dim(OUT$spectrum)
OUT$marginal
}
