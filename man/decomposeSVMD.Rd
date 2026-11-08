\name{decomposeSVMD}
\alias{decomposeSVMD}
\encoding{UTF-8}
\title{Successive Variational Mode Decomposition}
\description{
Extracts variational modes one at a time from the whole signal, keeping the
centers of the modes already extracted fixed, and stops after a known
number of modes or when the residual power meets a known noise power.
}
\usage{
decomposeSVMD(signal, samplingInterval, initialCenter = 0, alpha, maxAlpha,
              tau, tolerance, maxIterations = 299, stopping)
}
\arguments{
\item{signal}{Non-empty numeric vector of finite values without \code{dim},
in the units of the quantity.}
\item{samplingInterval}{Sampling interval \eqn{\Delta t} in seconds, a
finite number greater than zero.}
\item{initialCenter}{Initial center of every mode in Hz, in
\eqn{[0, 1/(2\Delta t)]}; default \code{0}.}
\item{alpha}{Initial penalty \eqn{\alpha > 0} in \eqn{s^2}. No default.}
\item{maxAlpha}{Final penalty, at least \code{alpha}, in \eqn{s^2}; equal to
\code{alpha} gives a fixed penalty. No default. PySDKit uses 10 and 20000
for its denominator \eqn{1 + 2 \alpha_\nu (\nu - \nu_L)^2} on cycles per
sample \eqn{\nu}; a value \eqn{\alpha_\nu} there corresponds to
\eqn{\alpha_\nu \Delta t^2 / (4 \pi^2)}.}
\item{tau}{Dual ascent step \eqn{\tau \ge 0}. No default; PySDKit uses 0.}
\item{tolerance}{Threshold \eqn{\epsilon_1 > 0} of the relative change of
the mode spectrum between inner iterations. No default; PySDKit uses
\code{1e-6}.}
\item{maxIterations}{Whole number \eqn{\ge 1} of inner iterations for each
value of the penalty; default \code{299}. The defaults of
\code{initialCenter} and \code{maxIterations} are those of the authors'
svmd.m 1.1.1.}
\item{stopping}{Named list with the field \code{stopRule} and exactly the
fields of that rule:
\describe{
\item{\code{"modeCount"}}{\code{modeCount}, the whole number of modes to
extract, without default.}
\item{\code{"noisePower"}}{\code{noisePower}, the noise power
\eqn{\sigma^2 > 0} in squared units of the signal; \code{noiseTolerance},
the dimensionless \eqn{\epsilon_2 > 0}; \code{maxModes}, a whole number
\eqn{\ge 1} or \code{Inf}, default \code{Inf}.}
}}
}
\details{
The signal is extended and transformed as in \code{\link{decomposeVMD}}.
Mode \eqn{L} is extracted from the whole signal with the centers
\eqn{\omega_1, \ldots, \omega_{L-1}}{omega_1, ..., omega_{L-1}} of the previous modes fixed. Its
penalty starts at \code{alpha} and doubles while below \code{maxAlpha};
when doubling would exceed \code{maxAlpha}, one last computation runs at
\code{maxAlpha}. Each penalty restarts the inner iteration from the
spectrum and center reached, with a zero multiplier and the penalty
\eqn{H_{L-1}(\omega) = \sum_{i<L} 1 / (\alpha^2 (\omega - \omega_i)^4)}. With
\eqn{A = \alpha^2 (\omega - \omega_L)^4}, the inner update is
\deqn{\hat u_L \leftarrow \frac{\hat x + A \hat u_L + \hat\lambda / 2}{(1 + A)\left(1 + 2 \alpha (\omega - \omega_L)^2 + H_{L-1}\right)},}{u_L <- (x + A u_L + lambda / 2) / ((1 + A) (1 + 2 alpha (omega - omega_L)^2 + H_{L-1})),}
followed by the new center and the dual ascent, and it stops when the
relative change of \eqn{\hat u_L}{u_L} is below \code{tolerance} or after
\code{maxIterations}. Extraction stops with \code{"modeCount"} after
\code{modeCount} modes; with \code{"noisePower"} when
\eqn{|\sigma^2 - P_L| / \sigma^2 < \epsilon_2}, \eqn{P_L} the mean square of
the residual on the record; or with \code{"modeLimit"} after
\code{maxModes} modes without meeting it. With \code{maxModes = Inf} the
extraction under \code{"noisePower"} continues until the condition holds.
The modes keep the extraction order. Their bandwidth formulation does not
impose the geometrical IMF conditions of EMD.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes, in extraction order.}
\item{fitResidual}{Numeric vector \eqn{x - \sum_L u_L}.}
\item{centers}{Numeric vector of the centers in Hz.}
\item{iterations}{Integer vector by mode, the inner iterations of all its
penalties.}
\item{converged}{Logical vector by mode, \code{TRUE} when every inner loop
of the mode met the change test.}
\item{stopReason}{\code{"modeCount"}, \code{"noisePower"} or
\code{"modeLimit"}.}
}
\section{Conditions}{
A mode without spectral energy leaves its center undefined and stops the
function with an error of class \code{xemdZeroEnergy}.
}
\references{
Nazari, M. and Sakhaei, S. M. (2020). Successive variational mode
decomposition. \emph{Signal Processing}, 107610.
\doi{10.1016/j.sigpro.2020.107610}
}
\seealso{
\code{\link{decomposeVMD}}, \code{\link{decomposeMVMD}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
OUT <- decomposeSVMD(Signal, samplingInterval = SamplingInterval,
                     initialCenter = 0,
                     alpha = 10 * SamplingInterval^2 / (4 * pi^2),
                     maxAlpha = 20000 * SamplingInterval^2 / (4 * pi^2),
                     tau = 0, tolerance = 1e-6, maxIterations = 299,
                     stopping = list(stopRule = "modeCount", modeCount = 5))
OUT$centers
OUT$stopReason
}
