\name{computeInstantaneousAttributes}
\alias{computeInstantaneousAttributes}
\encoding{UTF-8}
\title{Instantaneous Amplitude, Phase and Frequency}
\description{
Computes, for each column of a matrix of modes, the analytic signal and its
instantaneous amplitude, continuous phase and frequency.
}
\usage{
computeInstantaneousAttributes(modes, samplingInterval)
}
\arguments{
\item{modes}{Finite double matrix, samples by modes, with at least three
rows; it may have zero columns. Attributes other than \code{dim} do not
reach the outputs.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
}
\details{
The analytic signal \eqn{z = c + i \mathcal{H} c}{z = c + i H c} of each column is
computed with the discrete Fourier transform, multiplying each bin by
\eqn{1 + \mathrm{sgn}(\omega)}{1 + sgn(omega)}: positive frequencies are doubled, negative
ones removed, and the zero bin and, for an even number of samples, the
Nyquist bin keep the factor 1, so \eqn{\mathrm{Re}\, z = c}{Re z = c}. The transform
implies the periodic extension of the record. The amplitude is
\eqn{a = |z|}. Within each run of samples with \eqn{|z| > 0}, the phase
starts at the principal value of \eqn{\arg z}{arg z} and adds the increments
\eqn{\arg(z_{n+1} \bar z_n)}{arg(z_{n+1} conj(z_n))} taken in \eqn{(-\pi, \pi]}; a zero of \eqn{z}
leaves its sample undefined (\code{NA}) and the next run starts again at
its principal value. The frequency is
\deqn{f_n = \frac{\phi_{n+1} - \phi_{n-1}}{4 \pi \Delta t}}{f_n = (phi_{n+1} - phi_{n-1}) / (4 pi dt)}
in Hz, the central difference of \eqn{\omega = \phi'}; it is \code{NA} where
one of the three phases is, and each border repeats its neighbour. These
local undefined values raise no condition.
}
\value{
A list of three matrices, samples by modes:
\item{amplitude}{Instantaneous amplitude, in the units of the modes.}
\item{phase}{Continuous phase in radians, \code{NA} at zeros of the
analytic signal.}
\item{frequency}{Instantaneous frequency in Hz, \code{NA} where the phase
is undefined.}
}
\seealso{
\code{\link{buildHilbertSpectrum}}, \code{\link{decomposeEMD}},
\code{\link{estimateMaskFrequency}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
# Three tones of 2, 24 and 288 Hz as columns.
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
OUT <- computeInstantaneousAttributes(Tones,
                                      samplingInterval = SamplingInterval)
apply(OUT$amplitude, 2, range)
apply(OUT$frequency, 2, range)
}
