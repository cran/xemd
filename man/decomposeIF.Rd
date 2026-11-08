\name{decomposeIF}
\alias{decomposeIF}
\encoding{UTF-8}
\title{Iterative Filtering}
\description{
Decomposes a signal by successive extraction of components, each obtained by
iterating the subtraction of a moving average whose filter is a supplied
prototype resampled to a support set by the extrema of the residue.
}
\usage{
decomposeIF(signal, filter, boundary, tolerance, maxIterations, stretch,
            maxModes = Inf)
}
\arguments{
\item{signal}{Numeric vector of finite values without \code{dim}, in the
units of the quantity. Positions are sample indices.}
\item{filter}{Prototype filter, a numeric vector of finite weights of odd
length, at least 3. No default. The package tests use the 13 weights that
fifpy obtains from its prototype for a half width of 6.}
\item{boundary}{\code{"zero"} extends the record with zeros;
\code{"periodic"} extends it periodically and evaluates the iteration
spectrally (fast iterative filtering). No default.}
\item{tolerance}{Bound of the relative change between consecutive
iterates, greater than zero and less than one. No default; the package
tests use \eqn{\sqrt{0.001}}{sqrt(0.001)}.}
\item{maxIterations}{Whole number \eqn{\ge 1} of iterations of one
extraction. No default; the package tests use 200.}
\item{stretch}{Number greater than zero that scales the support of the
filter. No default; the package tests use 1.6.}
\item{maxModes}{Whole number \eqn{\ge 1} or \code{Inf}, the limit of
extractions; default \code{Inf}.}
}
\details{
Before extraction \eqn{k} the decomposition ends with
\code{"terminalResidue"} when \eqn{r_{k-1}} has at most one extremum and
with \code{"modeLimit"} when \eqn{k - 1} equals \code{maxModes}. Otherwise
the half width is \eqn{l_k = 2\,\mathrm{round}(\mathrm{stretch}\, N / E)}{l_k = 2 round(stretch N / E)}
for the \eqn{E} extrema of \eqn{r_{k-1}}; from the second extraction a value
that does not exceed the previous half width becomes
\eqn{\lceil 1.1\, l_{k-1} \rceil}{ceiling(1.1 l_{k-1})}. The prototype is resampled to
\eqn{2 l_k + 1} weights and the iteration
\deqn{f^{(q)} = f^{(q-1)} - \mathcal{L}_w f^{(q-1)}, \qquad f^{(0)} = r_{k-1},}{f^(q) = f^(q-1) - L_w f^(q-1), f^(0) = r_{k-1},}
with \eqn{\mathcal{L}_w}{L_w} the moving average under \code{boundary}, stops at
the first \eqn{q} with
\eqn{\|f^{(q)} - f^{(q-1)}\|_2 / \|f^{(q-1)}\|_2 \le}{||f^{(q)} - f^{(q-1)}||_2 / ||f^{(q-1)}||_2 <=} \code{tolerance}
(\code{converged} \code{TRUE}) or after \code{maxIterations}
(\code{converged} \code{FALSE}). Then \eqn{c_k = f^{(Q)}} and
\eqn{r_k = r_{k-1} - c_k}. A shrinking prototype is integrated over
rescaled bins; an enlarged one is interpolated linearly and normalized to a
unit absolute sum. The components are not shown to satisfy the counting
conditions of intrinsic mode functions.
}
\value{
A list with
\item{modes}{Double matrix, samples by components.}
\item{residue}{Numeric vector, the final residue.}
\item{iterations}{Integer vector, the iterations \eqn{Q} of each
component.}
\item{converged}{Logical vector by component.}
\item{stopReason}{\code{"terminalResidue"} or \code{"modeLimit"}.}
}
\references{
Cicone, A. and Zhou, H. (2018). Numerical Analysis for Iterative Filtering
with New Efficient Implementations Based on FFT. arXiv:1802.01359.

The support rule, the resampling of the prototype and both iterations are
adapted from fifpy, revision 26e51e7d45e551950dea8be37ece447fc4340296,
copyright (c) 2022 Emanuele Papini, MIT License; see the NOTICE file of the
installed package.
}
\seealso{
\code{\link{decomposeLMD}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq_len(1000) / 1000
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
Filter <- c(6.6306270449357867e-06, 0.00053827687799438354,
            0.0082858777882808092, 0.044245208647455882, 0.11927281750086954,
            0.20513332487239222, 0.24503568511284984, 0.20513332487239222,
            0.11927281750086952, 0.044245208647455903, 0.0082858777882808109,
            0.00053827687799438213, 6.6306270449358087e-06)
OUT <- decomposeIF(Signal, filter = Filter, boundary = "periodic",
                   tolerance = sqrt(0.001), maxIterations = 200,
                   stretch = 1.6, maxModes = 3)
OUT$iterations
OUT$stopReason
Zero <- decomposeIF(Signal, filter = Filter, boundary = "zero",
                    tolerance = sqrt(0.001), maxIterations = 200,
                    stretch = 1.6, maxModes = 3)
Zero$iterations
}
