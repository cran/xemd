\name{assessWhiteNoise}
\alias{assessWhiteNoise}
\encoding{UTF-8}
\title{White-Noise Significance Test of Modes}
\description{
Compares the energy of each mode with reference bounds from a white-noise
model at the mode's mean period (test of Wu and Huang).
}
\usage{
assessWhiteNoise(modes, probabilities, bound, rescalingMode)
}
\arguments{
\item{modes}{Finite numeric matrix, samples by modes. The caller supplies
the components of a normalized signal; see Details for the assumptions
of the reference distribution.}
\item{probabilities}{Two numbers \eqn{0 < p_- < p_+ < 1}, the
probabilities of the lower and upper bounds. No default; PyEMD uses
\eqn{\alpha = 0.95}, that is \eqn{(1 \mp \alpha)/2}{(1 -+ alpha)/2}.}
\item{bound}{\code{"gaussian"}, the line of the logarithmic quantile, or
\code{"chiSquare"}, the quantile of the energy. No default; PyEMD has only
the gaussian bound.}
\item{rescalingMode}{\code{NULL} for the test without rescaling, or the
index \eqn{r} of the column whose energy is placed on its upper bound. No
default; PyEMD rescales on the first component.}
}
\details{
The energy of mode \eqn{k} is \eqn{E_k = \sum_n c_k[n]^2 / N} and its mean
period \eqn{T_k} is \eqn{N} divided by its number of maxima, or \eqn{N}
without maxima. The bound of probability \eqn{p} is, with
\code{"gaussian"} and \eqn{z_p} the standard normal quantile,
\deqn{\ln E_p = -\ln T + z_p \sqrt{2/N}\, T^{1/2},}{ln E_p = -ln T + z_p sqrt(2 / N) T^(1/2),}
and with \code{"chiSquare"} \eqn{E_p = q(p; r) / N}, the chi-square
quantile with \eqn{r = N / T} degrees of freedom. With
\code{rescalingMode = r}, the energy of component \eqn{r} is placed on its
upper bound and every other energy is multiplied by
\eqn{E_{p_+}(T_r) / E_r}. A mode lies outside when its energy is below the
lower bound or above the upper one.

The reference distribution is conditional on the noise model,
normalization and EMD procedure. Interpreting the rescaled test assumes
that the component selected by \code{rescalingMode} is dominated by noise;
its mode index does not establish that assumption. Exceeding
a bound does not identify a physical source or prove that a component is
free of noise. A white-noise interpretation of an acceleration record does
not transfer unchanged to modes of its integrated displacement.
}
\value{
A list of vectors with one element per mode:
\item{energy}{The energies, rescaled when \code{rescalingMode} is given.}
\item{period}{The mean periods, in samples.}
\item{lower, upper}{The bounds at those periods.}
\item{outside}{Logical, \code{energy < lower | energy > upper}.}
}
\section{Conditions}{
A column of zero energy stops the function with an error of class
\code{xemdZeroEnergy}.
}
\seealso{
\code{\link{assessIMF}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
# The four summands of the signal, divided by its population deviation.
Summands <- cbind(sin(pi * Time), sin(2 * pi * Time), sin(6 * pi * Time),
                  0.5 * Time) / sqrt(mean((Signal - mean(Signal))^2))
assessWhiteNoise(Summands, probabilities = c(0.025, 0.975),
                 bound = "gaussian", rescalingMode = NULL)$outside
assessWhiteNoise(Summands, probabilities = c(0.025, 0.975),
                 bound = "chiSquare", rescalingMode = 1)$outside
}
