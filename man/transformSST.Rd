\name{transformSST}
\alias{transformSST}
\encoding{UTF-8}
\title{Synchrosqueezed Wavelet Transform}
\description{
Computes the Morlet continuous wavelet transform of a sampled signal and
reassigns its coefficients to frequency cells by the instantaneous
frequency estimated from the transform.
}
\usage{
transformSST(signal, samplingInterval, voiceCount = 16,
             morletCenter = 2 * pi, threshold)
}
\arguments{
\item{signal}{Numeric vector of at least four finite values without
\code{dim}, in the units of the quantity.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
\item{voiceCount}{Whole number \eqn{\ge 1} of voices per octave; default
\code{16}, the value of SynchWave.}
\item{morletCenter}{Center \eqn{\mu > 0} of the Morlet wavelet in
normalized frequency; default \eqn{2\pi}, the value of SynchWave.}
\item{threshold}{Number greater than zero in the units of the signal: only
coefficients with \eqn{|W| >} \code{threshold} are reassigned. No default;
SynchWave uses \code{sqrt(.Machine$double.eps)}.}
}
\details{
The record is extended by symmetric padding to
\eqn{T = 2^{1 + \mathrm{round}(\log_2 N)}}{T = 2^{1 + round(log2 N)}} samples, and the transform uses
the scales \eqn{a_j = 2^{j / \mathrm{voiceCount}}}{a_j = 2^{j / voiceCount}},
\eqn{j = 1, \ldots, (\log_2 T - 1)\,\mathrm{voiceCount}}{j = 1, ..., (log2 T - 1) voiceCount}. For each
coefficient above the threshold, the frequency
\deqn{\frac{\mathrm{Re}\,\Omega_f(a, b)}{2\pi} = \frac{1}{2\pi} \mathrm{Im}\frac{\partial_b W_f(a, b)}{W_f(a, b)}}{Re Omega / (2 pi) = Im(dW / W) / (2 pi)}
in Hz goes to the nearest of \eqn{J} logarithmic cells from
\eqn{1/((N - 1)\Delta t)} to \eqn{1/(2\Delta t)}, \eqn{J} the number of
scales, where the coefficient adds \eqn{W a^{-1/2}}; the sums are divided by
\code{voiceCount}. Estimates that are not positive or fall outside the
cells are not accumulated, and a transform whose energy is below the
machine epsilon gives a zero result. The representation recovers the
signal only approximately. The threshold is in the units of the signal,
so \code{\link{decomposeTSL}} has no operation for this transform.
}
\value{
A list with
\item{transform}{Complex matrix, samples by frequency cells.}
\item{frequency}{Numeric vector of the cell frequencies in Hz.}
}
\references{
Daubechies, I., Lu, J. and Wu, H.-T. (2009). Synchrosqueezed wavelet
transforms: a tool for empirical mode decomposition. arXiv:0912.2437.

The wavelet transform, the Morlet wavelet and the reassignment are adapted
from SynchWave 1.1.2, revision 9ac80b2c3e60b417d0ee5567a66e8e64faa1a979
(MATLAB original by Eugene Brevdo; R port by Dongik Jang, Hee-Seok Oh and
Donghoh Kim; the Synchrosqueezing toolbox attribution also names Hau-Tieng
Wu), LGPL-2.0-or-later, which \code{xemd} distributes under GPL-3 as section 3 of
the GNU Library General Public License version 2 permits; see the NOTICE
file of the installed package.
}
\seealso{
\code{\link{decomposeEWT}}, \code{\link{computeInstantaneousAttributes}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
OUT <- transformSST(Signal, samplingInterval = SamplingInterval,
                    voiceCount = 16, morletCenter = 2 * pi,
                    threshold = sqrt(.Machine$double.eps))
dim(OUT$transform)
range(OUT$frequency)
}
