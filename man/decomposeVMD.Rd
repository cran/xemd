\name{decomposeVMD}
\alias{decomposeVMD}
\encoding{UTF-8}
\title{Variational Mode Decomposition}
\description{
Decomposes a sampled signal into a given number of modes, each concentrated
around a center frequency, by alternating updates of the mode spectra, the
centers and a multiplier.
}
\usage{
decomposeVMD(signal, samplingInterval, initialCenters, alpha, tau,
             fixedZeroCenter, tolerance, maxIterations)
}
\arguments{
\item{signal}{Non-empty numeric vector of finite values without \code{dim},
in the units of the quantity.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
\item{initialCenters}{Numeric vector of initial center frequencies in Hz,
each in \eqn{[0, 1/(2\Delta t)]}; its length is the number of modes
\eqn{K}. No default.}
\item{alpha}{Penalty \eqn{\alpha > 0} of the bandwidth, in \eqn{s^2}
because it multiplies squared angular frequencies in rad/s. No default.
sktime and VMDecomp take \eqn{\alpha} on cycles per sample without the
factor 2: their value \eqn{\alpha_\nu} corresponds to
\eqn{\alpha = \alpha_\nu \Delta t^2 / (8 \pi^2)}.}
\item{tau}{Dual ascent step \eqn{\tau \ge 0}, dimensionless. \code{0}
solves the penalized problem with a zero multiplier; a positive value the
constrained problem. No default; the package tests use 0 and 0.5.}
\item{fixedZeroCenter}{\code{TRUE} keeps the first center at 0 and requires
\code{initialCenters[1] == 0}; \code{FALSE} updates every center. No
default.}
\item{tolerance}{Threshold \eqn{\epsilon > 0} of the relative change of the
mode spectra between sweeps, dimensionless. No default; the \code{1e-7} of
sktime belongs to an absolute criterion, not to this relative change.}
\item{maxIterations}{Whole number \eqn{\ge 1} of sweeps. No default;
sktime and VMDecomp allocate 500 states.}
}
\details{
The record of \eqn{N} samples is extended to \eqn{T = 2N} samples by
reflecting \eqn{\lceil N/2 \rceil}{ceiling(N/2)} samples at the start and the rest at the
end, and transformed with the unnormalized discrete Fourier transform; the
grid holds \eqn{|\omega|} in rad/s. Sweeps start from zero mode spectra and
multiplier. In each sweep mode \eqn{k}, updated after modes
\eqn{1, \ldots, k - 1}{1, ..., k - 1}, is
\deqn{\hat u_k = \frac{\hat x - \sum_{i \ne k} \hat u_i + \hat\lambda / 2}{1 + 2 \alpha (\omega - \omega_k)^2},}{u_k = (x - sum_{i != k} u_i + lambda / 2) / (1 + 2 alpha (omega - omega_k)^2),}
its center becomes the power-weighted mean of \eqn{|\omega|}, and the sweep
closes with \eqn{\hat\lambda \leftarrow \hat\lambda + \tau (\hat x - \sum_k \hat u_k)}{lambda <- lambda + tau (x - sum_k u_k)}.
The iteration stops with \code{"tolerance"} when
\eqn{\sum_k \|\hat u_k^{(n+1)} - \hat u_k^{(n)}\|^2 / \|\hat u_k^{(n)}\|^2 < \epsilon}{sum_k ||u_k^{(n+1)} - u_k^{(n)}||^2 / ||u_k^{(n)}||^2 < epsilon},
which cannot happen before the second sweep, or with
\code{"iterationLimit"} after \code{maxIterations} sweeps; the last sweep
is returned. The modes are the inverse transforms restricted to the record.
With \code{tau = 0} the penalized formulation permits a nonzero fit residual.
VMD defines its modes through a bandwidth objective, without imposing the
geometrical IMF conditions of EMD.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes, on the input grid.}
\item{fitResidual}{Numeric vector \eqn{x - \sum_k u_k}.}
\item{centers}{Numeric vector of the \eqn{K} centers of the last sweep, in
Hz.}
\item{iterations}{Integer vector repeating, for every mode, the number of
sweeps.}
\item{converged}{Logical vector repeating, for every mode, the result of
the joint stopping test.}
\item{stopReason}{\code{"tolerance"} or \code{"iterationLimit"}.}
}
\section{Conditions}{
A mode without spectral energy leaves its center undefined and stops the
function with an error of class \code{xemdZeroEnergy}.
}
\references{
Dragomiretskiy, K. and Zosso, D. (2014). Variational Mode Decomposition.
\emph{IEEE Transactions on Signal Processing}, 62(3), 531--544.
\doi{10.1109/TSP.2013.2288675}
}
\seealso{
\code{\link{decomposeMVMD}}, \code{\link{decomposeSVMD}},
\code{\link{decomposeTSL}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
OUT <- decomposeVMD(Signal, samplingInterval = SamplingInterval,
                    initialCenters = (0:2) / (6 * SamplingInterval),
                    alpha = 2000 * SamplingInterval^2 / (8 * pi^2),
                    tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6,
                    maxIterations = 98)
OUT$centers
OUT$stopReason
max(abs(OUT$fitResidual))
}
