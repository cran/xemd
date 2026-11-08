\name{decomposeMVMD}
\alias{decomposeMVMD}
\encoding{UTF-8}
\title{Multivariate Variational Mode Decomposition}
\description{
Decomposes the channels of a sampled multichannel signal into modes whose
channel coordinates share one center frequency.
}
\usage{
decomposeMVMD(signal, samplingInterval, initialCenters, alpha, tau,
              fixedZeroCenter, tolerance, maxIterations = 499)
}
\arguments{
\item{signal}{Numeric matrix of finite values, samples by channels, with at
least one row and one column, as composed by the caller and without
standardization.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
\item{initialCenters}{Numeric vector of initial center frequencies in Hz,
each in \eqn{[0, 1/(2\Delta t)]}; its length is the number of modes
\eqn{K}. No default.}
\item{alpha}{Penalty \eqn{\alpha > 0} in \eqn{s^2}, as in
\code{\link{decomposeVMD}}. No default; PySDKit takes it on cycles per
sample without the factor 2.}
\item{tau}{Dual ascent step \eqn{\tau \ge 0}; \code{0} solves the penalized
problem with zero multipliers. No default.}
\item{fixedZeroCenter}{\code{TRUE} keeps the first center at 0, while its
channel coordinates are still updated, and requires
\code{initialCenters[1] == 0}. No default.}
\item{tolerance}{Threshold \eqn{\epsilon > 0} of the relative change summed
over modes and channels. No default; the \code{1e-7} of PySDKit belongs to
another measure.}
\item{maxIterations}{Whole number \eqn{\ge 1} of sweeps; default
\code{499}, the limit of 500 states of the authors' MVMD.m 1.0.1.}
}
\details{
Each channel is extended and transformed as in \code{\link{decomposeVMD}}.
In each sweep the coordinate \eqn{c} of mode \eqn{k} is
\deqn{\hat u_{k,c} = \frac{\hat x_c - \sum_{i \ne k} \hat u_{i,c} + \hat\lambda_c / 2}{1 + 2 \alpha (\omega - \omega_k)^2},}{u_kc = (x_c - sum_{i != k} u_ic + lambda_c / 2) / (1 + 2 alpha (omega - omega_k)^2),}
with a denominator common to the channels; the new center weights
\eqn{|\omega|} by the power summed over the channels, and the multiplier of
each channel is updated with step \eqn{\tau}. The iteration stops with
\code{"tolerance"} when
\eqn{\sum_k \sum_c \|\hat u_{k,c}^{(n+1)} - \hat u_{k,c}^{(n)}\|^2 / \|\hat u_{k,c}^{(n)}\|^2 < \epsilon}{sum_k sum_c ||u_{k,c}^{(n+1)} - u_{k,c}^{(n)}||^2 / ||u_{k,c}^{(n)}||^2 < epsilon},
which cannot happen before the second sweep, or with
\code{"iterationLimit"} after \code{maxIterations} sweeps. A coordinate
without energy in the previous sweep leaves its ratio undefined, so the
test fails until \code{maxIterations}. The channels keep the order
received. The joint bandwidth formulation does not impose the geometrical
IMF conditions of EMD. A shared spectral center does not require identical
amplitudes or instantaneous frequencies across channels.
}
\value{
A list with
\item{modes}{Double array, samples by modes by channels.}
\item{fitResidual}{Numeric matrix, samples by channels,
\eqn{x_c - \sum_k u_{k,c}}.}
\item{centers}{Numeric vector of the \eqn{K} shared centers in Hz.}
\item{iterations}{Integer vector repeating, for every mode, the number of
sweeps.}
\item{converged}{Logical vector repeating, for every mode, the result of
the joint stopping test.}
\item{stopReason}{\code{"tolerance"} or \code{"iterationLimit"}.}
}
\section{Conditions}{
A mode without energy over all its channels leaves its center undefined and
stops the function with an error of class \code{xemdZeroEnergy}.
}
\references{
ur Rehman, N. and Aftab, H. (2019). Multivariate Variational Mode
Decomposition. arXiv:1907.04509.
}
\seealso{
\code{\link{decomposeVMD}}, \code{\link{decomposeSVMD}},
\code{\link{decomposeTSL}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cbind(cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
                    0.0625 * cos(2 * pi * 288 * Time),
                cos(2 * pi * 2 * Time),
                0.25 * cos(2 * pi * 24 * Time))
OUT <- decomposeMVMD(Signal, samplingInterval = SamplingInterval,
                     initialCenters = (0:2) / (6 * SamplingInterval),
                     alpha = 2000 * SamplingInterval^2 / (8 * pi^2),
                     tau = 0.5, fixedZeroCenter = FALSE, tolerance = 1e-6,
                     maxIterations = 98)
dim(OUT$modes)
OUT$centers
OUT$stopReason
}
