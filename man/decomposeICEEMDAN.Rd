\name{decomposeICEEMDAN}
\alias{decomposeICEEMDAN}
\encoding{UTF-8}
\title{Improved Complete Ensemble Empirical Mode Decomposition with Adaptive Noise}
\description{
Computes each residue as the average local mean of the previous residue
perturbed by a mode of supplied noise realizations, and each mode as the
difference between consecutive residues.
}
\usage{
decomposeICEEMDAN(signal, noise, amplitudeRatio = 0.2, normalizeNoise,
                  sifting)
}
\arguments{
\item{signal}{Numeric vector of at least two finite values without
\code{dim}, in the units of the quantity.}
\item{noise}{Finite numeric matrix with one row per sample of \code{signal}
and one column per realization, \eqn{M \ge 1}; column \eqn{m} is a
realization \eqn{w^{(m)}} of unit variance, used as supplied. The function
draws no random numbers.}
\item{amplitudeRatio}{Dimensionless ratio \eqn{\varepsilon_0 \ge 0}{epsilon_0 >= 0}.
Default \code{0.2}, as in \code{\link{decomposeEEMD}}.}
\item{normalizeNoise}{\code{TRUE} or \code{FALSE}, without default. From
the second stage on, \code{TRUE} divides each noise mode by its standard
deviation before scaling; \code{FALSE} keeps its amplitude.
\code{normalizeNoise = TRUE} corresponds to \code{SNRFlag = 2} of
ceemdan.m.}
\item{sifting}{Named list of sifting controls, as in
\code{\link{decomposeEMD}}; the rule \code{"consecutive"} is rejected. The
same list decomposes the noise realizations.}
}
\details{
Each column of \code{noise} is decomposed once with \code{sifting},
\eqn{\eta_k^{(m)} = E_k(w^{(m)})}. From \eqn{r_0 = x}, before stage \eqn{k}
the common residue \eqn{r_{k-1}} ends the decomposition, in this order, with
\code{"terminalResidue"} when it has at most one extremum, with
\code{"modeLimit"} when \eqn{k - 1} equals \code{sifting$maxModes} and with
\code{"insufficientExtrema"} when its envelopes cannot be built. Stage
\eqn{k} perturbs \eqn{r_{k-1}} with
\eqn{p^{(m)} = \beta_{k-1}^{(m)} \eta_k^{(m)}}, where
\eqn{\beta_0^{(m)} = \varepsilon_0 \, \mathrm{sd}(x) / \mathrm{sd}(\eta_1^{(m)})}{beta_0^{(m)} = epsilon_0 sd(x) / sd(eta_1^{(m)})}
and, for \eqn{k \ge 2}, \eqn{\beta_{k-1} = \varepsilon_0 \, \mathrm{sd}(r_{k-1})}{beta_{k-1} = epsilon_0 sd(r_{k-1})},
divided by \eqn{\mathrm{sd}(\eta_k^{(m)})}{sd(eta_k^{(m)})} when \code{normalizeNoise} is
\code{TRUE}; standard deviations use the divisor \eqn{N - 1}. With
\eqn{A} the average first mode of \eqn{r_{k-1} + p^{(m)}} and \eqn{\bar p}{mean(p)}
the mean perturbation, the new residue is the mean local mean and the mode
the difference of residues:
\deqn{r_k = r_{k-1} + \bar p - A, \qquad d_k = r_{k-1} - r_k.}{r_k = r_{k-1} + mean(p) - A, d_k = r_{k-1} - r_k.}
A realization without \eqn{\eta_k} contributes a null perturbation and
stays in the divisor \eqn{M}.
}
\value{
A list with
\item{modes}{Double matrix, samples by modes.}
\item{residue}{Numeric vector, the final common residue.}
\item{iterations}{Integer vector by mode, the sum of the subtractions of
the \eqn{M} extractions of the stage.}
\item{converged}{Logical vector by mode: \code{TRUE} when the \eqn{M}
extractions accepted, \code{FALSE} when one ended without accepting or
could not begin, \code{NA} when all completed fixed sifting.}
\item{stopReason}{\code{"terminalResidue"}, \code{"modeLimit"} or
\code{"insufficientExtrema"}.}
}
\section{Conditions}{
A noise mode of zero standard deviation leaves its scale undefined and
stops the function with an error of class \code{xemdZeroEnergy}.
}
\references{
Colominas, M. A., Schlotthauer, G. and Torres, M. E. (2014). Improved
complete ensemble EMD: A suitable tool for biomedical signal processing.
\emph{Biomedical Signal Processing and Control}, 14, 19--29.
\doi{10.1016/j.bspc.2014.06.009}

The selection of the reflected extrema and the spline interpolation of the
envelope nodes are adapted from emd.m (G. Rilling, March 2007), distributed
with CEEMDAN revision af9a0f83, copyright (c) 2020 Marcelo A. Colominas, MIT
License; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeCEEMDAN}}, \code{\link{decomposeEEMD}},
\code{\link{decomposeCEEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 1000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
set.seed(628)
Noise <- t(matrix(stats::rnorm(10L * length(Signal)), nrow = 10L,
                  byrow = TRUE))
Sifting <- list(stopRule = "doubleThreshold", thresholds = c(0.05, 0.5),
                fraction = 0.05, maxSift = 2000, boundary = "reflection",
                maxModes = Inf)
# Two of the ten realizations keep the example short.
OUT <- decomposeICEEMDAN(Signal, noise = Noise[, 1:2], amplitudeRatio = 0.2,
                         normalizeNoise = TRUE, sifting = Sifting)
dim(OUT$modes)
OUT$stopReason
max(abs(rowSums(OUT$modes) + OUT$residue - Signal))
}
