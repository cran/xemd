\name{decomposeMaskedEMD}
\alias{decomposeMaskedEMD}
\encoding{UTF-8}
\title{Masked Empirical Mode Decomposition}
\description{
Decomposes a signal by successive extraction in which each mode is the
average first mode of the residue plus sinusoidal masks of several phases,
under one of three mask rules: the opposite sine pair of 2005, dyadic masks
(EMD-DMS) and uniformly distributed phases (UPEMD).
}
\usage{
decomposeMaskedEMD(signal, samplingInterval, masking, sifting)
}
\arguments{
\item{signal}{Numeric vector of finite real values without \code{dim}, in
the units of the quantity.}
\item{samplingInterval}{Sampling interval in seconds, a finite number
greater than zero, without default. Mask frequencies are in Hz.}
\item{masking}{Named list with the field \code{rule} and exactly the fields
of that rule, without defaults:
\describe{
\item{\code{"pair"}}{\code{frequency}, in Hz, in
\eqn{(0, 1/(2\Delta t)]}, for instance the value of
\code{\link{estimateMaskFrequency}}; \code{amplitude}, greater than zero,
in the units of the signal.}
\item{\code{"dyadic"}}{\code{firstFrequency}, in Hz, in
\eqn{(0, 1/(2\Delta t)]}; \code{amplitudeRatio}, the dimensionless
\eqn{V_m > 0}. The example of the EMD-DMS article uses \eqn{V_m = 10}.}
\item{\code{"uniformPhase"}}{\code{amplitudeRatio}, the dimensionless
\eqn{\varepsilon_0 > 0}{epsilon_0 > 0}; \code{phaseCount}, the whole number \eqn{P \ge 2}
of phases. The package tests use \eqn{\varepsilon_0 = 1}{epsilon_0 = 1} and
\eqn{P = 4}.}
}}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected.}
}
\details{
From \eqn{r_0 = x}, before stage \eqn{k} the decomposition ends, in this
order, with \code{"terminalResidue"} when \eqn{r_{k-1}} has at most one
extremum, with \code{"modeLimit"} when \eqn{k - 1} equals
\code{sifting$maxModes} or the schedule of masks is exhausted, and with
\code{"insufficientExtrema"} when the envelopes of \eqn{r_{k-1}} cannot be
built. Otherwise the masks of stage \eqn{k} are
\eqn{a_k \cos(2 \pi f_k n \Delta t + \theta_p)}{a_k cos(2 pi f_k n Delta t + theta_p)}, \eqn{n = 0, \ldots, N - 1}{n = 0, ..., N - 1},
with the first sample as phase origin, and the mode is the average of the
first sifted component of \eqn{r_{k-1}} plus each mask,
\deqn{c_k = \frac{1}{P} \sum_{p=1}^{P} E_1\left(r_{k-1} + w_p\right), \qquad r_k = r_{k-1} - c_k.}{c_k = (1/P) sum_p E_1(r_{k-1} + w_p), r_k = r_{k-1} - c_k.}
A stage is retained unless one of its phases has no first mode or loses
its envelopes during sifting; then the decomposition stops before adding
the stage mean or subtracting it from the residue. A phase that reaches
the iteration limit is kept and flags the stage \code{FALSE}. The phases and divisor are unchanged.
The rules fix the phases, frequencies and amplitudes:
\describe{
\item{\code{"pair"}}{Phases \eqn{-\pi/2} and \eqn{\pi/2}, so the masks
are \eqn{\pm a \sin(2 \pi f n \Delta t)}{+- a sin(2 pi f n Delta t)} and the mode is
\eqn{(E_1(x + w) + E_1(x - w)) / 2}. The schedule has one stage: the
decomposition ends with \code{"modeLimit"} after one extraction unless the
residue is terminal.}
\item{\code{"dyadic"}}{Phases 0 and \eqn{\pi};
\eqn{f_k = 2^{1-k} f_1}; amplitude \eqn{V_m} times the mean absolute
value of the maxima of the residue.}
\item{\code{"uniformPhase"}}{Phases \eqn{2 \pi p / P},
\eqn{p = 0, \ldots, P - 1}{p = 0, ..., P - 1}; \eqn{f_k = 1 / (2^k \Delta t)}; amplitude
\eqn{\varepsilon_0}{epsilon_0} times the standard deviation of the residue with
divisor \eqn{N}.}
}
The masks are deterministic: the function draws no random numbers.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes.}
\item{residue}{Numeric vector, the final residue.}
\item{iterations}{Integer vector by mode, the sum of the subtractions of
the extractions of the stage.}
\item{converged}{Logical vector by mode: \code{TRUE} when every extraction
of the stage accepted; \code{NA} when all completed fixed sifting.
Interrupted stages are reported in \code{termination}.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"},
\code{"insufficientExtrema"} or \code{"incompleteExtraction"}, the last
when at least one phase fails.}
\item{maskFrequencies}{Numeric vector, the mask frequency of each mode in
Hz.}
\item{termination}{Data frame with \code{mode}, the next stage index after
the retained modes; \code{phase}, its one-based index in the phase sequence;
\code{iterations}, subtractions of that failed phase; and \code{stopReason},
\code{"terminalResidue"}, \code{"modeLimit"}, \code{"insufficientExtrema"}
or \code{"iterationLimit"}. One row is returned for each phase of the interrupted stage that did not
accept, including the phases that reached the iteration limit.
If the decomposition stops before a stage starts, one row carries its
\code{stopReason}, \code{phase = NA} and zero iterations.}
}
\references{
Deering, R. and Kaiser, J. F. (2005). The use of a masking signal to
improve empirical mode decomposition. In \emph{IEEE International
Conference on Acoustics, Speech, and Signal Processing (ICASSP)}, 4,
485--488.

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{estimateMaskFrequency}},
\code{\link{decomposeIteratedMaskedEMD}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = 3)
# Uniform phases, with one sample per second.
OUT <- decomposeMaskedEMD(Signal, samplingInterval = 1,
                          masking = list(rule = "uniformPhase",
                                         amplitudeRatio = 1, phaseCount = 4),
                          sifting = Sifting)
OUT$maskFrequencies
OUT$stopReason

# The pair of 2005 at the estimated mask frequency.
Frequency <- estimateMaskFrequency(Signal, samplingInterval = 1,
                                   sifting = Sifting)
Pair <- decomposeMaskedEMD(Signal, samplingInterval = 1,
                           masking = list(rule = "pair", frequency = Frequency,
                                          amplitude = sqrt(mean((Signal -
                                              mean(Signal))^2))),
                           sifting = Sifting)
dim(Pair$modes)
Pair$stopReason
}
