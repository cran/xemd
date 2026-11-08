\name{decomposeTVFEMD}
\alias{decomposeTVFEMD}
\encoding{UTF-8}
\title{Time-Varying Filter Empirical Mode Decomposition}
\description{
Decomposes a signal by successive extraction of components, each obtained by
sifting with a local mean that is a B-spline approximation on knots placed
by a time-varying cutoff frequency (TVF-EMD).
}
\usage{
decomposeTVFEMD(signal, bandwidthThreshold, splineOrder,
                intermittencyThreshold, maxSift, maxModes = Inf)
}
\arguments{
\item{signal}{Numeric vector of finite real values without \code{dim}, in
the units of the quantity. Positions are sample indices and the cutoff is
in radians per sample.}
\item{bandwidthThreshold}{Dimensionless threshold \eqn{\xi > 0} of the
relative instantaneous bandwidth. No default; the original study uses
0.1.}
\item{splineOrder}{Whole number \eqn{n \ge 1}, the order of the B-spline
(degree \eqn{n - 1}). No default; the original study uses 28.}
\item{intermittencyThreshold}{Dimensionless threshold \eqn{\rho > 0} of
the intermittency test of the cutoff. No default; the original study uses
0.25.}
\item{maxSift}{Whole number \eqn{\ge 1} of subtractions of one component.
No default; pytvfemd stops after five.}
\item{maxModes}{Whole number \eqn{\ge 1} or \code{Inf}, the limit of
components; default \code{Inf}.}
}
\details{
Before extraction \eqn{k} the decomposition ends with
\code{"terminalResidue"} when \eqn{r_{k-1}} has fewer than four extrema,
maxima plus minima, and with \code{"modeLimit"} when \eqn{k - 1} equals
\code{maxModes}. A component is obtained by sifting,
\deqn{h_{s+1} = h_s - m_s,}{h_{s+1} = h_s - m_s,}
where \eqn{m_s} is the periodic least-squares B-spline approximation of
order \eqn{n} of \eqn{h_s} on knots at the extrema of
\eqn{\cos\varphi}{cos phi}, \eqn{\varphi}{phi} the cumulative realigned cutoff. Each
iteration extends \eqn{h_s} by symmetric reflection of
\eqn{\mathrm{round}(N/2)}{round(N/2)} samples per side and restores the record
afterwards. The cutoff follows from the amplitude \eqn{A} and the frequency
of the analytic signal, interpolated at the minima and at the maxima of
\eqn{A}; values above
\eqn{0.9\pi} are lowered to \eqn{0.9\pi}, undefined or non-positive values
are interpolated, and the zones marked by intermittency events, where the
relative variation of the cutoff between consecutive maxima of the signal
exceeds \eqn{\rho}, are replaced by interpolation.

The stopping rule aggregates the relative instantaneous bandwidth
\eqn{\theta} over the samples from \eqn{\mathrm{round}(0.05N)}{round(0.05N)} to
\eqn{\mathrm{round}(0.95N)}{round(0.95N)}: \eqn{h_s} is accepted when the mean of
\eqn{\theta} is below \eqn{\xi (1 + (s + 1)/4)}, from \eqn{s = 1} in the
first component and from \eqn{s = 0} in the later ones, or when the samples
with \eqn{\theta > \xi} are fewer than \eqn{0.2N}. This aggregation of
pytvfemd 0.3.14 replaces the local criterion of the monograph.
\code{converged} is \code{FALSE} when \code{maxSift} subtractions end
without acceptance or when the cutoff or the knots cannot be built.
Acceptance refers to the aggregated rule above; it does not establish the
pointwise narrow-band condition. Neither that rule nor a returned component
guarantees the two classical IMF conditions.
}
\value{
A list with
\item{modes}{Double matrix, samples by components.}
\item{residue}{Numeric vector, the final residue.}
\item{iterations}{Integer vector, the subtractions of each component.}
\item{converged}{Logical vector by component.}
\item{stopReason}{\code{"terminalResidue"} or \code{"modeLimit"}.}
}
\references{
The stopping rule and the reflection of the sifting, the end of the
decomposition, the interpolation, the bounds and the realignment of the
cutoff, the relative bandwidth, the knots and the periodic B-spline
approximation are adapted from pytvfemd 0.3.14 (Stefano Bianchi and
Alessandro Longo), GPL-3; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeEMD}}, \code{\link{computeInstantaneousAttributes}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
OUT <- decomposeTVFEMD(Signal, bandwidthThreshold = 0.1, splineOrder = 28,
                       intermittencyThreshold = 0.25, maxSift = 5)
dim(OUT$modes)
OUT$iterations
OUT$stopReason
}
