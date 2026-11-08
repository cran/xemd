\name{assessFractionalNoise}
\alias{assessFractionalNoise}
\encoding{UTF-8}
\title{Fractional Gaussian Noise Test of Intrinsic Mode Functions}
\description{
Compares the energy of each intrinsic mode function, from the second on,
with the threshold that fractional Gaussian noise of a tabulated Hurst
exponent gives at a tabulated level, with the first function as the noise
scale.
}
\usage{
assessFractionalNoise(modes, hurst, probability)
}
\arguments{
\item{modes}{Finite numeric matrix of 2048 samples by at least two
intrinsic mode functions, in their order of extraction and without the
residue.}
\item{hurst}{Hurst exponent \eqn{H}: \code{0.2}, \code{0.5} or \code{0.8}.
No default.}
\item{probability}{Level \eqn{p} of the threshold: \code{0.95} or
\code{0.99}. No default.}
}
\details{
The energy of function \eqn{k} is \eqn{W[k] = \sum_n c_k[n]^2}. From
\eqn{k = 2} on, the reference energy is
\deqn{W_H[k] = \frac{W[1]}{\beta_H} \rho_H^{-2(1-H)k}, \quad \rho_H = 2.01 + 0.2 (H - 0.5) + 0.12 (H - 0.5)^2,}{W_H[k] = W[1] / beta_H rho_H^(-2 (1 - H) k), rho_H = 2.01 + 0.2 (H - 0.5) + 0.12 (H - 0.5)^2,}
and the threshold is
\eqn{T_{H,p}[k] = W_H[k] \, 2^{2^{a_{H,p} k + b_{H,p}}}}{T_{H,p}[k] = W_H[k] 2^(2^(a_{H,p} k + b_{H,p}))}, with
\eqn{\beta_H}, \eqn{a_{H,p}} and \eqn{b_{H,p}} of table 3.1 of Flandrin et
al. (2005). A function lies outside when \eqn{W[k] > T_{H,p}[k]}. The
coefficients come from simulations of 10000 realizations of 2048 samples
and are not interpolated to other exponents nor recalibrated for other
lengths. The first function is taken to carry mainly noise; the function
does not verify it.

The result is conditional on that noise model and calibration. Exceeding
the threshold does not prove that a component is free of noise or identify
its physical source.
}
\value{
A list of vectors with one element per function:
\item{energy}{The energies \eqn{W[k]}.}
\item{reference}{The reference energies \eqn{W_H[k]}, \code{NA} for the
first function.}
\item{upper}{The thresholds \eqn{T_{H,p}[k]}, \code{NA} for the first
function.}
\item{outside}{Logical, \code{energy > upper}, \code{NA} for the first
function.}
}
\section{Conditions}{
A first function of zero energy stops the function with an error of class
\code{xemdZeroEnergy}.
}
\references{
Flandrin, P., \enc{Gonçalvès}{Goncalves}, P. and Rilling, G. (2004).
Detrending and denoising with empirical mode decompositions.
\emph{EUSIPCO 2004}, 1581--1584.

Flandrin, P., \enc{Gonçalvès}{Goncalves}, P. and Rilling, G. (2005). EMD
equivalent filter banks, from interpretation to applications. In N. E.
Huang and S. S. P. Shen (eds.), \emph{Hilbert-Huang Transform and Its
Applications}, World Scientific, 57--74.
}
\seealso{
\code{\link{assessWhiteNoise}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq_len(2048) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
assessFractionalNoise(Tones, hurst = 0.5, probability = 0.95)
# The third tone as the noise scale.
assessFractionalNoise(Tones[, c(3, 1, 2)], hurst = 0.5,
                      probability = 0.95)$outside
}
