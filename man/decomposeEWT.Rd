\name{decomposeEWT}
\alias{decomposeEWT}
\encoding{UTF-8}
\title{Empirical Wavelet Transform}
\description{
Splits the spectrum of a sampled signal at the midpoints between its
dominant spectral maxima, or at the boundaries of a detection rule of the
EWT code of Gilles, and returns the components of an adapted filter bank,
which add up to the signal.
}
\usage{
decomposeEWT(signal, samplingInterval, bandCount, detection = NULL)
}
\arguments{
\item{signal}{Numeric vector of at least four finite values without
\code{dim}, in the units of the quantity.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
\item{bandCount}{Whole number \eqn{\ge 2}, the largest number of bands, or
\code{NULL} for the automatic count of \code{detection}. No default; ewtpy
uses 5, and pyewt, the package of the author of the method, requires the
value.}
\item{detection}{\code{NULL}, the default, for the midpoint boundaries, or
a named list with these fields:
\describe{
\item{\code{detector}}{\code{"locmax"}, \code{"locmaxmin"} or
\code{"locmaxminf"}; required.}
\item{\code{logarithm}}{\code{TRUE} to work on the natural logarithm of the
spectral magnitude, which must then be positive; default \code{FALSE}.}
\item{\code{regularization}}{\code{"none"}, \code{"average"} or
\code{"gaussian"}; default \code{"gaussian"}.}
\item{\code{window}}{Whole number \eqn{L \ge 1} of cells of the
regularization filter; default \code{10}.}
\item{\code{sigma}}{Standard deviation \eqn{\sigma > 0} of the gaussian
filter, in cells; default \code{1.5}.}
\item{\code{thresholdRatio}}{Number \eqn{\alpha} in \eqn{(0, 1)}, given
only with \code{bandCount = NULL}. No default; Gilles (2013) tries 0.3 and
0.4.}
}
The defaults are the initial values of the EWT code of Gilles; its initial
detector, \code{"scalespace"}, is not available. With \code{bandCount = NULL} the list
must give \code{detector = "locmax"}, \code{regularization = "none"},
\code{logarithm = FALSE} and \code{thresholdRatio}.}
}
\details{
The maxima of the one-sided magnitude of the discrete Fourier transform are
ranked, and the \code{bandCount - 1} largest are retained, or all of them
when fewer exist, so the number of bands adapts to the spectrum. With
\eqn{b_0 = 0} and \eqn{b_1 < \cdots < b_{K-1}}{b_1 < ... < b_{K-1}} the bins of the retained
maxima, the boundaries are the midpoints
\eqn{\Omega_j = \pi (b_{j-1} + b_j) / N} rad per sample, and the last band
closes at \eqn{\pi}. The record is extended by mirroring
\eqn{\lceil N/2 \rceil - 1}{ceiling(N/2) - 1} samples before it and \eqn{\lceil N/2 \rceil}{ceiling(N/2)}
after it. On that extension the bank holds an approximation filter and
band filters with cosine transitions of half width \eqn{\gamma \Omega_j},
\eqn{\gamma} below the bound that keeps neighbouring transitions apart. The
analysis coefficients are
\eqn{W_j = \mathcal{F}^{-1}[\hat x \, \overline{H_j}]}{W_j = IDFT(DFT(x) conj(H_j))}, and each component
applies the dual filter \eqn{H_j / \sum_i H_i^2},
\deqn{c_j = \mathcal{F}^{-1}\left[\frac{\hat W_j H_j}{\sum_i H_i^2}\right],}{c_j = IDFT(DFT(W_j) H_j / sum_i H_i^2),}
restricted to the record, so the components add up to the signal. The
components are defined by the filter bank, without imposing the IMF
conditions of EMD. The method has no iterative stopping rule.

A detection list works on the one-sided magnitude \eqn{s}, or on its
natural logarithm, and on its regularization \eqn{p}, the central \eqn{n}
values of the full convolution of \eqn{s} with the weights \eqn{h_r},
\eqn{r = 0, \ldots, L-1}{r = 0, ..., L - 1}, from position
\eqn{\lfloor L/2 \rfloor + 1}{floor(L/2) + 1}: \eqn{h_r = 1/L} for
\code{"average"} and
\eqn{h_r \propto \exp[-u_r^2 / (2\sigma^2)]}{h_r proportional to exp(-u_r^2 / (2 sigma^2))},
\eqn{u_r = r - (L-1)/2}, for \code{"gaussian"}; \code{"none"} keeps
\eqn{p = s}. \code{"locmax"} places the midpoints above on the maxima of
\eqn{p}. \code{"locmaxmin"} and \code{"locmaxminf"} need
\code{bandCount - 1} positive maxima of \eqn{p} and retain the largest; on
the interval from the first cell to the first retained maximum and on those
between consecutive retained maxima, endpoints included, the boundary is
the lowest minimum of \eqn{p} for \code{"locmaxmin"} or of \eqn{s} for
\code{"locmaxminf"}, the upper central one among equal minima, and a
minimum at position \eqn{q} gives \eqn{\Omega = 2\pi (q - 1) / N}. No
boundary is completed. With \code{bandCount = NULL}, the number of bands is
the number of maxima of the magnitude whose height exceeds
\eqn{M_m + \alpha (M_1 - M_m)}{M_m + alpha (M_1 - M_m)}, \eqn{M_1} and
\eqn{M_m} the largest and smallest heights, and the midpoint boundaries
follow; one band returns the signal.
}
\value{
A list with
\item{modes}{Double matrix, samples by bands.}
\item{boundaries}{Numeric vector of the \eqn{K - 1} boundaries in Hz.}
}
\section{Conditions}{
Equal magnitudes at the last retained maximum and the first discarded one
leave the selection undefined and stop with an error, as does a frequency
where every filter of the bank vanishes. With a detection list, a zero
magnitude under \code{logarithm = TRUE}, fewer than \code{bandCount - 1}
positive maxima of \eqn{p} or an interval without minimum for
\code{"locmaxmin"} and \code{"locmaxminf"}, boundaries that do not increase
and an automatic count of zero bands stop with an error.
}
\references{
Gilles, J. (2013). Empirical wavelet transform. \emph{IEEE Transactions on
Signal Processing}, 61(16), 3999--4010.

Gilles, J. \emph{Empirical Wavelets}, MATLAB toolbox, revision
0381208d40352237876b79912013b0f68337133c,
\url{https://github.com/jegilles/Empirical-Wavelets}.

The detection of the boundaries, the mirror extension and the filter bank
are adapted from ewtpy, revision f9ff38bc5cd7859eea58c5b29099f8a503792b63,
copyright (c) 2019 \enc{Vinícius}{Vinicius} Carvalho and Eduardo Mazoni, MIT
License; the synthesis through the dual bank is adapted from pyewt, revision
9d60a340f34e92356ed4f7a6a8b7044f524455e3, copyright (c) 2024 Jerome Gilles,
MIT License. See the NOTICE file of the installed package.
}
\seealso{
\code{\link{transformSST}}, \code{\link{decomposeVMD}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
OUT <- decomposeEWT(Signal, samplingInterval = SamplingInterval,
                    bandCount = 3)
OUT$boundaries
dim(OUT$modes)
max(abs(rowSums(OUT$modes) - Signal))
# Automatic count: only the 2 Hz maximum exceeds the threshold of 0.3.
dim(decomposeEWT(Signal, samplingInterval = SamplingInterval,
                 bandCount = NULL,
                 detection = list(detector = "locmax", regularization = "none",
                                  thresholdRatio = 0.3))$modes)
# Lowest minima between the maxima of the regularized magnitude.
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
decomposeEWT(Signal, samplingInterval = 9 / 2999, bandCount = 3,
             detection = list(detector = "locmaxmin"))$boundaries
}
