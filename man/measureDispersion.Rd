\name{measureDispersion}
\alias{measureDispersion}
\encoding{UTF-8}
\title{Dispersion of Ensemble Members}
\description{
Computes, sample by sample and mode by mode, the sample variance of the
members of an ensemble around their mean, and that sample variance divided
by the number of members. The latter estimates the variance of the mean
under the assumptions stated below.
}
\usage{
measureDispersion(realizations)
}
\arguments{
\item{realizations}{Finite numeric array with three dimensions, samples by
modes by realizations, with at least two realizations. The modes of each
realization occupy their extraction order, and a mode absent from a
realization is filled with zeros.}
}
\details{
With \eqn{M} realizations \eqn{C^{(m)}} and their mean \eqn{\bar C}{mean(C)}, the
sample variance is
\deqn{s^2 = \frac{1}{M - 1} \sum_{m=1}^{M} \left(C^{(m)} - \bar C\right)^2}{s^2 = sum_m (C^(m) - mean(C))^2 / (M - 1)}
and the returned estimate for the variance of the mean is \eqn{s^2 / M}.
The first describes dispersion between members. Interpreting the second
as Monte Carlo precision assumes comparable components and independent
realizations with a common mean and a common finite covariance.
Opposite perturbations in a complementary pair are dependent; the two
signs do not provide independent realizations for this interpretation.

The function does not establish these statistical assumptions or physical
correspondence between equal mode indices. Dispersion under perturbations
of a fixed record does not quantify all uncertainty in the physical signal,
and increasing the number of members does not remove decomposition bias.
Attributes of \code{realizations} other than \code{dim} do not reach the
outputs.
}
\value{
A list with two matrices, samples by modes:
\item{variance}{The sample variance \eqn{s^2}, in squared units of the
signal.}
\item{meanVariance}{The estimate \eqn{s^2 / M}, subject to the assumptions
in Details.}
}
\seealso{
\code{\link{decomposeEEMD}}
}
\examples{
# Two realizations of four samples: the first has two modes, the second
# one, so its second mode is filled with zeros.
Realizations <- array(c(2, 4, 6, 8, 1, 1, 1, 1,
                        4, 6, 8, 10, 0, 0, 0, 0), dim = c(4, 2, 2))
measureDispersion(Realizations)
}
