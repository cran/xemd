\name{measureOrthogonality}
\alias{measureOrthogonality}
\encoding{UTF-8}
\title{Inner Products and Cosines Between Components}
\description{
Computes the inner products between the columns of a matrix of components
and the cosines of the angles between them.
}
\usage{
measureOrthogonality(components)
}
\arguments{
\item{components}{Finite numeric matrix, samples by components. Attributes
other than \code{dim} do not reach the outputs.}
}
\details{
The inner products are the Euclidean products \eqn{\sum_n v_i[n] v_j[n]},
without the sampling interval of the discrete inner product, which cancels
in the cosine
\deqn{\cos\theta_{ij} = \frac{\langle v_i, v_j \rangle}{\|v_i\|_2 \|v_j\|_2}.}{cos theta_ij = <v_i, v_j> / (||v_i|| ||v_j||).}
The cosine keeps its sign. The cross products are the terms by which the
energy of a sum of components differs from the sum of their energies.
}
\value{
A list of two square matrices, components by components:
\item{innerProduct}{The inner products.}
\item{cosine}{The cosines.}
}
\section{Conditions}{
A column of zero norm stops the function with an error of class
\code{xemdZeroEnergy}.
}
\seealso{
\code{\link{measurePMSI}}, \code{\link{measureFCF}}
}
\examples{
Time <- seq_len(1000) / 1000
Tones <- cbind(cos(2 * pi * 2 * Time), 0.25 * cos(2 * pi * 24 * Time),
               0.0625 * cos(2 * pi * 288 * Time))
# The sum of the tones and the tones.
OUT <- measureOrthogonality(cbind(rowSums(Tones), Tones))
round(OUT$cosine, 6)
}
