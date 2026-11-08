\name{measurePMSI}
\alias{measurePMSI}
\encoding{UTF-8}
\title{Pseudo-Mode Splitting Index}
\description{
Computes the pseudo-mode splitting index of each pair of adjacent modes.
}
\usage{
measurePMSI(modes, denominator)
}
\arguments{
\item{modes}{Finite numeric matrix, samples by modes, with at least one
column. Attributes other than \code{dim} do not reach the outputs.}
\item{denominator}{\code{"homogeneous"}, the sum of the squared norms of
the pair, or \code{"regularized"}, which adds the \eqn{10^{-8}} of the
UPEMD definition, an absolute scale absent from the homogeneous form. No
default; emd 0.8.1 has only the regularized form.}
}
\details{
For each pair of adjacent columns,
\deqn{\mathrm{PMSI}(c_i, c_{i+1}) = \max\left\{\frac{\langle c_i, c_{i+1} \rangle}{\|c_i\|_2^2 + \|c_{i+1}\|_2^2}, 0\right\},}{PMSI(c_i, c_{i+1}) = max(<c_i, c_{i+1}> / (||c_i||^2 + ||c_{i+1}||^2), 0),}
with \eqn{10^{-8}} added to the denominator under \code{"regularized"}.
}
\value{
Numeric vector with one value per pair of adjacent modes, of length one
less than the number of columns.
}
\section{Conditions}{
A pair of adjacent modes with zero energy under \code{"homogeneous"} stops
the function with an error of class \code{xemdZeroEnergy}.
}
\seealso{
\code{\link{measureOrthogonality}}, \code{\link{decomposeMaskedEMD}}
}
\examples{
Time <- seq_len(1000) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
Mixed <- cbind(rowSums(Tones), Tones)
measurePMSI(Mixed, denominator = "homogeneous")
measurePMSI(Mixed, denominator = "regularized")
}
