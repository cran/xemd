\name{decomposeMEMD}
\alias{decomposeMEMD}
\encoding{UTF-8}
\title{Multivariate Empirical Mode Decomposition}
\description{
Decomposes a signal of two or more channels by successive extraction of
vector modes, with vector envelopes built at the maxima of its projections
onto directions of the channel space.
}
\usage{
decomposeMEMD(signal, directionCount = 64, sifting)
}
\arguments{
\item{signal}{Numeric matrix of finite real values, samples by channels,
with two or more columns, as composed by the caller and without
standardization.}
\item{directionCount}{Whole number \eqn{D \ge 1} of Hammersley points;
default \code{64}, common to Rlibeemd and emdr. The projections use the
\eqn{2D} unit vectors formed by those points and their antipodes.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}, with the rules \code{"fixed"},
\code{"consecutive"} and \code{"doubleThreshold"} and
\code{boundary = "reflection"} only. A missing \code{stopRule} is
\code{"doubleThreshold"}, whose \code{thresholds}, \code{fraction} and
\code{maxSift} take the defaults of \code{\link{decomposeEMD}}; its extra
field \code{amplitudeRule}, \code{"globalMean"} or \code{"envelopePairs"}
(default), selects the amplitude of the ratio. Under \code{"consecutive"},
\code{consecutiveCount} (a count, default 2) and \code{maxSift} (default
1000, at least \code{consecutiveCount}) apply. The default rule, the
default amplitude and the defaults of \code{"consecutive"} are those of
memd.m by Rehman and Mandic.}
}
\details{
The directions map \eqn{D} Hammersley points, built from the radical
inverses in the bases of the first \eqn{C - 1} primes, to unit vectors
\eqn{\mathbf{v}_d}{v_d} through the angles of \code{get.projections} of emdr
1.2.3 (Pierre Masselot); their antipodes complete \eqn{2D} directions. For each
direction, \eqn{\mathbf{e}_d}{e_d} interpolates the vector samples at the maxima
of the projection \eqn{p_d = \mathbf{v}_d^{T} \mathbf{x}}{p_d = v_d^{T} x}, extended by
reflection, and the sifting subtracts the mean of the \eqn{2D} envelopes,
\deqn{\mathbf{h}_{s+1} = \mathbf{h}_s - \frac{1}{2D} \sum_{d=1}^{2D} \mathbf{e}_d[\mathbf{h}_s].}{h_{s+1} = h_s - (1/(2D)) sum_d e_d[h_s].}
The adaptive rules act on the projections. \code{"consecutive"} accepts when
every projection keeps \eqn{|E - Z| \le 1} in \code{consecutiveCount}
successive states after a subtraction, the run restarting when one does
not. \code{"doubleThreshold"} accepts a state, the signal included, whose
projections all keep \eqn{|E - Z| \le 1} and whose ratio
\eqn{|\mathbf{m}| / a}{|m| / a} meets the double threshold of
\code{\link{decomposeEMD}}; \eqn{a} is the mean distance of the envelopes
to their mean with \code{"globalMean"}, or the mean half distance between
the envelopes of antipodal directions with \code{"envelopePairs"}. Before
extraction \eqn{k} the decomposition ends with \code{"terminalResidue"}
when every projection of \eqn{\mathbf{r}_{k-1}}{r_{k-1}} has at most one extremum,
with \code{"modeLimit"} when \eqn{k - 1} equals \code{sifting$maxModes}, and
with \code{"insufficientExtrema"} when a direction has no envelope. One
extraction serves all the channels, so they share the number of modes.
}
\value{
A list with
\item{modes}{Double array, samples by modes by channels, in the order of
the columns of \code{signal}.}
\item{residue}{Numeric matrix, samples by channels.}
\item{iterations}{Integer vector, the subtractions of each mode.}
\item{converged}{Logical vector by mode, as in \code{\link{decomposeEMD}}.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}.}
}
\references{
Rehman, N. and Mandic, D. P. (2010). Multivariate empirical mode
decomposition. \emph{Proceedings of the Royal Society A}, 466,
1291--1302. \doi{10.1098/rspa.2009.0502}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeNAMEMD}}, \code{\link{decomposeBEMD}},
\code{\link{decomposeMVMD}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- 2 * pi * (0:255) / 256
z <- cos(0.3 * Time) * exp(2i * Time) +
    0.3 * abs(sin(2.3 * Time)) * exp(17i * Time)
Signal <- cbind(Re(z), Im(z))
OUT <- decomposeMEMD(Signal, directionCount = 64,
                     sifting = list(stopRule = "fixed", fixedSift = 10,
                                    boundary = "reflection", maxModes = 3))
dim(OUT$modes)
OUT$stopReason
}
