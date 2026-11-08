\name{selectModesEnergy}
\alias{selectModesEnergy}
\encoding{UTF-8}
\title{Selection of Modes by Energy}
\description{
Orders the modes by decreasing energy and selects those before the largest
drop between consecutive energies.
}
\usage{
selectModesEnergy(modes)
}
\arguments{
\item{modes}{Finite numeric matrix, samples by modes, with at least two
columns. Attributes other than \code{dim} do not reach the outputs.}
}
\details{
The energy of mode \eqn{k} is \eqn{P_k = \sum_n c_k[n]^2}. With the
permutation \eqn{\pi} that orders the energies decreasingly, the selection
keeps \eqn{\pi(1), \ldots, \pi(v)}{pi(1), ..., pi(v)} with
\deqn{v = \arg\max_{1 \le m < K} \left(P_{\pi(m)} - P_{\pi(m+1)}\right).}{v = argmax over m of P_pi(m) - P_pi(m+1).}
The source fixes no tie-break, so a largest drop that is zero or not
unique stops the function with an error.
}
\value{
Integer vector of the selected column indices, the most energetic first.
}
\seealso{
\code{\link{measureFCF}}, \code{\link{measureDTW}}
}
\examples{
Time <- seq_len(1000) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
selectModesEnergy(Tones)
# The same tones in another order.
selectModesEnergy(Tones[, c(3, 1, 2)])
}
