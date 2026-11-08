\name{decomposeBEMD}
\alias{decomposeBEMD}
\encoding{UTF-8}
\title{Bivariate Empirical Mode Decomposition}
\description{
Decomposes a two-channel signal, read as the complex signal
\eqn{z = x_1 + i x_2}, by successive extraction of bivariate modes with
envelopes built on projections onto equally spaced directions of the plane.
}
\usage{
decomposeBEMD(signal, directionCount = 64, envelope, sifting)
}
\arguments{
\item{signal}{Numeric matrix of finite real values with two columns, the
real and the imaginary parts, samples in rows.}
\item{directionCount}{Whole number \eqn{D \ge 1} of directions
\eqn{\varphi_d = 2 \pi (d - 1) / D}{phi_d = 2 pi (d - 1) / D}; default \code{64}, the default of
Rlibeemd \code{bemd} and of emdr \code{memd}.}
\item{envelope}{Mean of the sifting, without default:
\code{"trajectory"} averages the envelopes of the trajectory itself at the
maxima of each projection; \code{"projection"} sums the envelopes of the
projections, oriented on their directions, with the factor \eqn{2/D}.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}, with \code{stopRule = "fixed"} only. The
trajectory envelope admits only \code{boundary = "reflection"}; the
projection envelope also admits \code{"extrapolation"}.}
}
\details{
The projection on direction \eqn{\varphi_d}{phi_d} is
\eqn{p_d = x_1 \cos\varphi_d + x_2 \sin\varphi_d}{p_d = x_1 cos phi_d + x_2 sin phi_d}. With
\code{"trajectory"}, \eqn{e_d} interpolates the two-dimensional samples at
the maxima of \eqn{p_d} and the mean is
\deqn{m = \frac{1}{D} \sum_{d=1}^{D} e_d;}{m = (1/D) sum_d e_d;}
with \code{"projection"}, the scalar upper envelope of \eqn{p_d} is
oriented on \eqn{(\cos\varphi_d, \sin\varphi_d)}{(cos phi_d, sin phi_d)} and the sum is scaled by
\eqn{2/D}. Each extraction subtracts the mean \code{fixedSift} times. Before
extraction \eqn{k} the decomposition ends with \code{"terminalResidue"}
when every projection of \eqn{r_{k-1}} has at most one extremum, with
\code{"modeLimit"} when \eqn{k - 1} equals \code{sifting$maxModes}, and with
\code{"insufficientExtrema"} when a direction has no envelope. One
extraction serves both channels, so they share the number of modes.
If a direction loses its envelope before fixed sifting completes, the
decomposition stops before adding that vector candidate and preserves
the residue from before the attempt. A candidate that reaches the
iteration limit of an adaptive rule is kept and flagged \code{FALSE}. Directions are not dropped from the mean.
}
\value{
A list with
\item{modes}{Double array, samples by modes by channels.}
\item{residue}{Numeric matrix, samples by channels.}
\item{iterations}{Integer vector, the subtractions of each mode.}
\item{converged}{Logical vector by mode: \code{NA} after the fixed number
of subtractions. Interrupted attempts are reported in \code{termination}.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}.}
\item{termination}{Data frame with one row: \code{mode}, the next extraction
index after the retained modes; \code{direction}, the one-based index of
the first direction without an envelope, or \code{NA} when none failed;
\code{iterations}, subtractions of the interrupted attempt, zero if it did
not start; and \code{stopReason}.}
}
\references{
Rilling, G., Flandrin, P., \enc{Gonçalves}{Goncalves}, P. and Lilly, J. M.
(2007). Bivariate Empirical Mode Decomposition. \emph{IEEE Signal
Processing Letters}, 14, 936--939.

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeMEMD}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- 2 * pi * (0:511) / 512
z <- cos(0.3 * Time) * exp(2i * Time) +
    0.3 * abs(sin(2.3 * Time)) * exp(17i * Time)
Signal <- cbind(Re(z), Im(z))
OUT <- decomposeBEMD(Signal, directionCount = 64, envelope = "projection",
                     sifting = list(stopRule = "fixed", fixedSift = 10,
                                    boundary = "extrapolation",
                                    maxModes = 3))
dim(OUT$modes)
OUT$stopReason
max(abs(apply(OUT$modes, c(1, 3), sum) + OUT$residue - Signal))
}
