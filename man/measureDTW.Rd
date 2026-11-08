\name{measureDTW}
\alias{measureDTW}
\encoding{UTF-8}
\title{Dynamic Time Warping Distance of Modes}
\description{
Computes the dynamic time warping distance between a signal and each of its
modes.
}
\usage{
measureDTW(signal, modes)
}
\arguments{
\item{signal}{Non-empty finite numeric vector without \code{dim}.}
\item{modes}{Finite numeric matrix with one row per sample of
\code{signal} and one column per mode. Attributes other than \code{dim} do
not reach the outputs.}
}
\details{
With the local cost \eqn{d_{ij} = (x_i - c_j)^2} and the accumulation
\deqn{G_{ij} = d_{ij} + \min\{G_{i-1,j-1}, G_{i-1,j}, G_{i,j-1}\},}{G_ij = d_ij + min(G_{i-1,j-1}, G_{i-1,j}, G_{i,j-1}),}
\eqn{G_{00} = 0} and infinite cost outside the matrix, the distance is
\eqn{\sqrt{G_{NN}}}{sqrt(G_{NN})}, without normalization by the length of the path. The
work grows with the square of the number of samples for each mode.
}
\value{
Numeric vector with one distance per mode, in the units of the signal.
}
\seealso{
\code{\link{measureFCF}}, \code{\link{selectModesEnergy}}
}
\examples{
Time <- seq_len(1000) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
measureDTW(signal = rowSums(Tones), modes = Tones)
}
