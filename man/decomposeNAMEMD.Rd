\name{decomposeNAMEMD}
\alias{decomposeNAMEMD}
\encoding{UTF-8}
\title{Noise-Assisted Multivariate Empirical Mode Decomposition}
\description{
Decomposes a signal together with supplied noise channels by one
multivariate decomposition and keeps the coordinates of the signal.
}
\usage{
decomposeNAMEMD(signal, noiseChannels, directionCount = 64, sifting)
}
\arguments{
\item{signal}{Numeric matrix of finite real values, samples by channels,
with one or more columns.}
\item{noiseChannels}{Numeric matrix of finite real values with the rows of
\code{signal} and \eqn{Q \ge 1} columns, in the units of the signal,
supplied by the caller: the function draws no random numbers. emdr draws
its noise with standard deviation 0.02 on standardized channels.}
\item{directionCount}{Whole number of Hammersley points, as in
\code{\link{decomposeMEMD}}; default \code{64}.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeMEMD}}.}
}
\details{
The function runs \code{\link{decomposeMEMD}} on the augmented matrix
\eqn{\mathbf{y} = [\mathbf{x}, \mathbf{w}]}{y = [x, w]} and keeps the first \eqn{C}
coordinates of every mode, \eqn{\mathbf{c}_k = P_x \tilde{\mathbf{c}}_k}{c_k = P_x c~_k},
and of the residue. Projections, stopping, iterations and termination are
those of the augmented decomposition. In \code{\link{decomposeTSL}} the
noise channels are drawn by group with \code{noiseChannelCount} and scaled
by \code{noiseRatio} times the reference deviation of the group.
}
\value{
A list with the fields of \code{\link{decomposeMEMD}}; \code{modes} has
the \eqn{C} channels of \code{signal} and \code{residue} its \eqn{C}
columns.
}
\references{
ur Rehman, N. and Mandic, D. P. (2011). Filter Bank Property of
Multivariate Empirical Mode Decomposition. \emph{IEEE Transactions on
Signal Processing}, 59(5), 2421--2426. \doi{10.1109/TSP.2011.2106779}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeMEMD}}, \code{\link{decomposeTSL}}
}
\examples{
Time <- 2 * pi * (0:511) / 512
z <- cos(0.3 * Time) * exp(2i * Time) +
    0.3 * abs(sin(2.3 * Time)) * exp(17i * Time)
Signal <- cbind(Re(z), Im(z))
set.seed(628)
# Two noise channels with 0.02 times the reference deviation of the signal.
Noise <- 0.02 * sqrt(mean(apply(Signal, 2, stats::sd)^2)) *
    matrix(stats::rnorm(512 * 2), nrow = 512)
OUT <- decomposeNAMEMD(Signal, noiseChannels = Noise, directionCount = 64,
                       sifting = list(stopRule = "fixed", fixedSift = 1,
                                      boundary = "reflection",
                                      maxModes = 3))
dim(OUT$modes)
dim(OUT$residue)
}
