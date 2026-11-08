\name{decomposeEMD}
\alias{decomposeEMD}
\encoding{UTF-8}
\title{Empirical Mode Decomposition}
\description{
Decomposes a signal by successive extraction: each component is obtained
by sifting the current residue and is subtracted from that residue.
Its numerical acceptance depends on the selected stopping rule.
}
\usage{
decomposeEMD(signal, sifting)
}
\arguments{
\item{signal}{Numeric vector of finite real values without \code{dim}, in
the units of the quantity. Positions are sample indices, so no time vector
or sampling interval is needed. Attributes such as names, \code{tsp} or
class do not reach the outputs.}
\item{sifting}{Named list of sifting controls: \code{stopRule} and the
fields of that rule, described in section \sQuote{Sifting controls}. Fields
with a default may be omitted; every other field is required, and a field
that does not belong to the rule is an error.}
}
\details{
With \eqn{r_0 = x}, extraction \eqn{k} computes the mode
\eqn{c_k = E_1(r_{k-1})} by sifting and updates the residue
\deqn{r_k = r_{k-1} - c_k,}{r_k = r_{k-1} - c_k,}
so that \eqn{x = \sum_{k=1}^{K} c_k + r_K} up to rounding. Before
extraction \eqn{k} the decomposition ends, in this order, with
\code{"terminalResidue"} when \eqn{r_{k-1}} has at most one extremum, with
\code{"modeLimit"} when \eqn{k - 1} equals \code{sifting$maxModes}, and with
\code{"insufficientExtrema"} when the boundary rule cannot build the
envelopes of \eqn{r_{k-1}}.

Sifting starts from \eqn{h_0 = r_{k-1}} and repeats
\eqn{h_{j+1} = h_j - m_j}, where \eqn{m_j} is the mean of the upper and
lower envelopes of \eqn{h_j}. The envelopes interpolate the maxima and the
minima, extended at both borders by \code{sifting$boundary}: a cubic spline
with not-a-knot end conditions through four or more nodes, the parabola
through three and the line through two. A run of equal interior samples
is one extremum at the exact centre of the run. Each state \eqn{h_j} is
examined once, in this order: the stopping rule accepts it
(\code{converged} \code{TRUE}); \eqn{j} equals \code{fixedSift} under
\code{"fixed"} (\code{converged} \code{NA}); \eqn{j} equals \code{maxSift}
under an adaptive rule, or the envelopes of \eqn{h_j} cannot be built
(\code{converged} \code{FALSE}); otherwise one more subtraction follows.
Non-convergence is reported in \code{converged}, without a warning.
}
\section{Sifting controls}{
The same list configures the extractor of \code{\link{decomposeEEMD}},
\code{\link{decomposeCEEMD}}, \code{\link{decomposeCEEMDAN}},
\code{\link{decomposeICEEMDAN}}, \code{\link{resiftModes}},
\code{\link{decomposeMaskedEMD}}, \code{\link{decomposeIteratedMaskedEMD}},
\code{\link{estimateMaskFrequency}} and, with the restrictions stated in
their help pages, \code{\link{decomposeBEMD}}, \code{\link{decomposeMEMD}}
and \code{\link{decomposeNAMEMD}}. Counts are accepted as doubles with a
whole value.
\describe{
\item{\code{stopRule}}{Categorical, without default: \code{"SD"},
\code{"Snumber"}, \code{"doubleThreshold"} or \code{"fixed"}. The rule
\code{"consecutive"} belongs to the vector sifting of
\code{\link{decomposeMEMD}} and is rejected here.}
\item{\code{tolerance}}{Under \code{"SD"}: dimensionless threshold,
greater than zero, without default. From \eqn{j = 1} the state is accepted
when \eqn{\sum_n (h_{j-1}[n] - h_j[n])^2 / h_{j-1}[n]^2} is finite and
below it; the sum runs over all samples, without the factor \eqn{1/N}, and
a null denominator is not regularized. The 1998 formulation proposes
values from 0.2 to 0.3.}
\item{\code{stableSift}}{Under \code{"Snumber"}: the count \eqn{S \ge 1},
without default. From \eqn{j = 1} the state is accepted when the counts of
extrema \eqn{E} and zero crossings \eqn{Z} have stayed equal over \eqn{S}
consecutive states and \eqn{|E - Z| \le 1}. The 2003 formulation reports
typical values from three to five.}
\item{\code{thresholds}}{Under \code{"doubleThreshold"}: two dimensionless
numbers \eqn{0 < \theta_1 < \theta_2}; default \code{c(0.05, 0.5)}.}
\item{\code{fraction}}{Under \code{"doubleThreshold"}: the admitted
fraction \eqn{\alpha} in \eqn{[0, 1)}; default \code{0.05}. With
\eqn{\sigma = |m / a|}, \eqn{m} the half sum and \eqn{a} the half
difference of the envelopes, the state, the signal included, is accepted
when \eqn{|E - Z| \le 1}, the fraction of samples with
\eqn{\sigma \ge \theta_1} is at most \eqn{\alpha} and no sample has
\eqn{\sigma \ge \theta_2}. Where \eqn{m = a = 0} the ratio is undefined and
counts as no exceedance; where only \eqn{a = 0} it is infinite.}
\item{\code{fixedSift}}{Under \code{"fixed"}: the count \eqn{J \ge 1} of
subtractions of every extraction, without default.}
\item{\code{maxSift}}{Under \code{"SD"}, \code{"Snumber"} and
\code{"doubleThreshold"}: the count \eqn{M \ge 1} of subtractions after
which an extraction ends without acceptance; under \code{"Snumber"}
\eqn{M \ge S}. Default \code{2000} under \code{"doubleThreshold"}, none
under the other rules. The rule \code{"fixed"} has no \code{maxSift}.}
\item{\code{maxModes}}{Count \eqn{\ge 1} or \code{Inf}, the limit of
extractions; default \code{Inf}, which extracts until the residue ends the
decomposition.}
\item{\code{boundary}}{Categorical: \code{"reflection"} or
\code{"extrapolation"}; default \code{"reflection"} under
\code{"doubleThreshold"} and \code{"fixed"}, none under \code{"SD"} and
\code{"Snumber"}. \code{"reflection"} mirrors two extrema of each family
per side about the first extremum or the border sample and needs three
extrema; \code{"extrapolation"} anchors the upper envelope at each border
at the larger of the border sample and the line through the two nearest
maxima, and the lower envelope at the smaller of the border sample and the
line through the two nearest minima.}
}
}
\value{
A list with
\item{modes}{Double matrix, samples by modes, in the units of
\code{signal}; it has zero columns when no mode is extracted. No names.}
\item{residue}{Numeric vector \eqn{r_K} of the input length; equal to
\code{signal} when no mode is extracted.}
\item{iterations}{Integer vector with the number of subtractions of each
mode.}
\item{converged}{Logical vector by mode: \code{TRUE} when the stopping rule
accepted the mode, \code{FALSE} when \code{maxSift} was reached or the
envelopes could not be built, \code{NA} under \code{"fixed"}.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}.}
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. C., Shih, H. H., Zheng, Q.,
Yen, N.-C., Tung, C. C. and Liu, H. H. (1998). The empirical mode
decomposition and the Hilbert spectrum for nonlinear and non-stationary
time series analysis. \emph{Proceedings of the Royal Society of London.
Series A}, 454, 903--995. \doi{10.1098/rspa.1998.0193}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{assessIMF}}, \code{\link{decomposeEEMD}},
\code{\link{decomposeTSL}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = Inf)
OUT <- decomposeEMD(Signal, sifting = Sifting)
dim(OUT$modes)
OUT$iterations
OUT$stopReason
max(abs(rowSums(OUT$modes) + OUT$residue - Signal))

# The rule alone: the other fields take their defaults.
identical(decomposeEMD(Signal, sifting = list(stopRule = "doubleThreshold")),
          OUT)

# Ten fixed subtractions per mode, at most three modes.
Fixed <- decomposeEMD(Signal, sifting = list(stopRule = "fixed",
                                             fixedSift = 10, maxModes = 3))
Fixed$converged
Fixed$stopReason
}
