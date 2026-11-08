\name{assessIMF}
\alias{assessIMF}
\encoding{UTF-8}
\title{Intrinsic Mode Function Conditions}
\description{
Evaluates numerical diagnostics of the intrinsic mode function conditions
for each column of a matrix: the balance of extrema and zero crossings,
and a double threshold on the envelope mean relative to its amplitude.
}
\usage{
assessIMF(modes, thresholds = c(0.05, 0.5), fraction = 0.05,
          boundary = "reflection")
}
\arguments{
\item{modes}{Finite numeric matrix, samples by columns. Attributes other
than \code{dim} do not reach the outputs.}
\item{thresholds}{Two dimensionless numbers \eqn{0 < \theta_1 < \theta_2};
default \code{c(0.05, 0.5)}.}
\item{fraction}{The admitted fraction \eqn{\alpha} in \eqn{[0, 1)};
default \code{0.05}.}
\item{boundary}{\code{"reflection"} (default) or \code{"extrapolation"},
the extension of the envelopes at the borders, as in
\code{\link{decomposeEMD}}. The defaults are those of emd.m, shared with
the sifting controls.}
}
\details{
The ideal envelope condition is a zero local mean. This diagnostic uses
tolerances and a specified boundary rule; passing it does not establish
recovery of a physical component. For methods defined by other conditions,
such as a bandwidth objective, it is an additional EMD-based diagnostic.

\code{count} tests \eqn{|N_e - N_z| \le 1}, with the extrema \eqn{N_e} and
zero crossings \eqn{N_z} counted as in \code{\link{decomposeEMD}}.
\code{deviation} applies to
\deqn{\sigma = \frac{|m_e|}{a_e},}{sigma = |m_e| / a_e,}
\eqn{m_e} and \eqn{a_e} the half sum and the half difference of the
envelopes, the double threshold that accepts a mode under
\code{stopRule = "doubleThreshold"}: the fraction of samples with
\eqn{\sigma \ge \theta_1} is at most \eqn{\alpha} and no sample has
\eqn{\sigma \ge \theta_2}.
}
\value{
A list of two logical vectors with one element per column:
\item{count}{The counting condition.}
\item{deviation}{The double threshold, \code{NA} where the boundary rule
leaves the column without envelopes.}
}
\references{
The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeEMD}}, \code{\link{assessWhiteNoise}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
OUT <- decomposeEMD(Signal, sifting = list(stopRule = "doubleThreshold"))
# The signal, its modes and its residue as columns.
assessIMF(cbind(Signal, OUT$modes, OUT$residue))
}
