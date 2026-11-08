\name{measureQRF}
\alias{measureQRF}
\encoding{UTF-8}
\title{Quality of Reconstruction Factor}
\description{
Computes the quality of reconstruction factor of an estimate against a
reference component, in decibels.
}
\usage{
measureQRF(estimate, reference)
}
\arguments{
\item{estimate}{Finite numeric vector without \code{dim}.}
\item{reference}{Finite numeric vector without \code{dim}, of the length of
\code{estimate}, with a nonzero norm.}
}
\details{
\deqn{\mathrm{QRF}(\hat s, s) = 20 \log_{10} \frac{\|s\|_2}{\|s - \hat s\|_2},}{QRF = 20 log10(||s|| / ||s - estimate||),}
\code{Inf} when the estimate equals the reference.
}
\value{
A number in dB.
}
\section{Conditions}{
A reference of zero norm is outside the domain of the factor and stops the
function with an error.
}
\seealso{
\code{\link{measureFCF}}, \code{\link{measureOrthogonality}}
}
\examples{
Time <- seq_len(1000) / 1000
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
# The signal without its third tone, against the signal.
measureQRF(estimate = Signal - 0.0625 * cos(2 * pi * 288 * Time),
           reference = Signal)
measureQRF(estimate = Signal, reference = Signal)
}
