\name{measureFCF}
\alias{measureFCF}
\encoding{UTF-8}
\title{Fault Correlation Factor of Modes}
\description{
Computes the fault correlation factor between a signal and each of its
modes, a normalized inner product without removal of means.
}
\usage{
measureFCF(signal, modes)
}
\arguments{
\item{signal}{Finite numeric vector without \code{dim}.}
\item{modes}{Finite numeric matrix with one row per sample of
\code{signal} and one column per mode. Attributes other than \code{dim} do
not reach the outputs.}
}
\details{
\deqn{\mathrm{FCF}(x, c) = \frac{\sum_n x[n] c[n]}{\sqrt{\sum_n x[n]^2 \sum_n c[n]^2}}.}{FCF(x, c) = sum(x c) / sqrt(sum(x^2) sum(c^2)).}
The means are not removed, so the factor differs from a correlation of
centered vectors when the inputs have a nonzero mean.
}
\value{
Numeric vector with one factor per mode.
}
\section{Conditions}{
A signal or a mode of zero norm stops the function with an error of class
\code{xemdZeroEnergy}.
}
\seealso{
\code{\link{measureDTW}}, \code{\link{selectModesEnergy}}
}
\examples{
Time <- seq_len(1000) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
measureFCF(signal = rowSums(Tones), modes = Tones)
}
