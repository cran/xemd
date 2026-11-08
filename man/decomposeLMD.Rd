\name{decomposeLMD}
\alias{decomposeLMD}
\encoding{UTF-8}
\title{Local Mean Decomposition}
\description{
Decomposes a signal by successive extraction of product functions, each
obtained by repeated subtraction of a smoothed local mean and division by a
smoothed local magnitude.
}
\usage{
decomposeLMD(signal, maxIterations, envelopeTolerance, changeTolerance,
             maxSmoothPasses, maxModes = Inf)
}
\arguments{
\item{signal}{Numeric vector of finite values without \code{dim}, in the
units of the quantity. Positions are sample indices.}
\item{maxIterations}{Whole number \eqn{\ge 1} of normalizations allowed in
one product function. No default; the package tests use 200.}
\item{envelopeTolerance}{Number \eqn{\ge 0}; an extraction ends when the
mean of \eqn{|1 - a_q|} over the samples is at most this value. No default;
the package tests use 0.01.}
\item{changeTolerance}{Number \eqn{\ge 0}; an extraction ends when the mean
of \eqn{|s_{q-1} - s_q|} is at most this value. No default; the package
tests use 0.01.}
\item{maxSmoothPasses}{Whole number \eqn{\ge 1} of passes of the moving
average that smooths the local levels. No default; 12 is the number of
passes that LMD 1.2.1 reaches by default.}
\item{maxModes}{Whole number \eqn{\ge 1} or \code{Inf}, the limit of
product functions; default \code{Inf}.}
}
\details{
The local levels of a function use as nodes its interior extrema and its
first and last samples: consecutive nodes \eqn{n_i, n_{i+1}} give the mean
\eqn{m_i = (n_i + n_{i+1}) / 2} and the magnitude
\eqn{a_i = |n_i - n_{i+1}| / 2}, held between the two nodes. Both levels are
smoothed by a repeated uniform moving average whose half width follows from
the largest gap between nodes and shrinks at the borders; passes repeat
while two consecutive samples are equal, up to \code{maxSmoothPasses}.

From \eqn{s_0 = r_{p-1}}, iteration \eqn{q} forms
\deqn{s_q = \frac{s_{q-1} - m_q}{a_q}}{s_q = (s_{q-1} - m_q) / a_q}
with the smoothed levels of \eqn{s_{q-1}}. The extraction ends without that
normalization, with \code{converged} \code{TRUE}, when one of the two
tolerances is met; it ends with \code{converged} \code{FALSE} when
\eqn{s_{q-1}} has at most one extremum or after \code{maxIterations}
normalizations. With \eqn{Q} normalizations, the product function is
\eqn{a_1 \cdots a_Q s_Q}{a_1 ... a_Q s_Q} and \eqn{r_p = r_{p-1} - PF_p}. Before extraction
\eqn{p} the decomposition ends with \code{"terminalResidue"} when
\eqn{r_{p-1}} has at most one extremum and with \code{"modeLimit"} when
\eqn{p - 1} equals \code{maxModes}. A product function combines an
accumulated amplitude envelope with a normalized oscillation; this
construction does not impose the IMF envelope conditions of EMD.

In the first iteration \eqn{a_1} and \eqn{s_0 = r_{p-1}} carry the units of
the signal, so both tolerances compare magnitudes in those units with
dimensionless quantities; the result depends on the unit only when the
initial magnitude approaches 1 in that unit.
}
\value{
A list with
\item{modes}{Double matrix, samples by product functions.}
\item{residue}{Numeric vector, the final residue.}
\item{iterations}{Integer vector, the number \eqn{Q} of normalizations of
each product function.}
\item{converged}{Logical vector by product function.}
\item{stopReason}{\code{"terminalResidue"} or \code{"modeLimit"}.}
}
\section{Conditions}{
A local magnitude that is not positive leaves the normalization undefined
and stops the function with an error of class \code{xemdZeroEnvelope}.
}
\references{
Smith, J. S. (2005). The local mean decomposition and its application to
EEG perception data. \emph{Journal of the Royal Society Interface}, 2(5),
443--454. \doi{10.1098/rsif.2005.0058}

The local levels and the extraction of a product function are adapted from
LMD 1.2.1 (Shubhra Prakash), revision
ec11b721ec7b05f0137eeb67a619b28f62ef85a2, Apache License (>= 2); see the
NOTICE file of the installed package.
}
\seealso{
\code{\link{decomposeIF}}, \code{\link{decomposeEMD}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
OUT <- decomposeLMD(Signal, maxIterations = 200, envelopeTolerance = 0.01,
                    changeTolerance = 0.01, maxSmoothPasses = 12,
                    maxModes = 8)
dim(OUT$modes)
OUT$iterations
OUT$stopReason
}
