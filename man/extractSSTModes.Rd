\name{extractSSTModes}
\alias{extractSSTModes}
\encoding{UTF-8}
\title{Components Along the Ridges of a Synchrosqueezed Transform}
\description{
Extracts ridges of a synchrosqueezed transform one after another by dynamic
programming and reconstructs a component from the coefficients in a band
around each ridge.
}
\usage{
extractSSTModes(transform, frequency, modeCount, penalty = 1000,
                exclusionWidth = 4, reconstructionWidth, morletCenter)
}
\arguments{
\item{transform}{Numeric or complex matrix of finite values, samples by
frequency cells, such as the \code{transform} of \code{\link{transformSST}}.}
\item{frequency}{Increasing positive frequencies of the cells, one per
column of \code{transform}, such as the \code{frequency} of
\code{\link{transformSST}}.}
\item{modeCount}{Whole number \eqn{\ge 1} of ridges. No default;
curve_ext_multi of the Synchrosqueezing Toolbox requires the value.}
\item{penalty}{Dimensionless penalty \eqn{\lambda \ge 0} of a jump between
cells; default \code{1000}, the value that curve_ext_multi assigns.}
\item{exclusionWidth}{Whole number \eqn{b \ge 0} of cells on each side of a
ridge that are removed before the next extraction; default \code{4}, the
value that curve_ext_multi assigns.}
\item{reconstructionWidth}{Whole number \eqn{w \ge 0} of cells on each side
of a ridge in its reconstruction band. No default; a value not below the
number of cells gives the complete inversion, as synsq_cwt_iw does without
a band.}
\item{morletCenter}{Center \eqn{\mu > 0} of the Morlet wavelet that
produced \code{transform}, in normalized frequency. No default.}
}
\details{
With \eqn{T_{n,j}} the coefficients, \eqn{f_j} the frequencies and
\eqn{z_j = \log_2 f_j}{z_j = log2(f_j)}, the cost of a cell is
\eqn{E_{n,j} = -\log(|T_{n,j}|^2 / Q + 10^{-8})}{E_{n,j} = -log(|T_{n,j}|^2 / Q + 1e-8)}, \eqn{Q} the sum of
\eqn{|T|^2} over the matrix. The accumulated cost is
\eqn{F_{1,j} = E_{1,j}} and
\deqn{F_{n,j} = E_{n,j} + \min_k \left[F_{n-1,k} + \lambda (z_j - z_k)^2\right]}{F_{n,j} = E_{n,j} + min_k [F_{n-1,k} + lambda (z_j - z_k)^2]}
over every cell \eqn{k}; the ridge ends at the first cell of least
\eqn{F_{N,j}}, and each earlier cell is the first whose residual in that
recurrence is below \eqn{10^{-8}} in modulus. After each ridge \eqn{C_n},
the cells from \eqn{C_n - b} to \eqn{C_n + b}, clipped to the cells, take
the value \eqn{\sqrt{\epsilon}}{sqrt(eps)}, \eqn{\epsilon} the machine
precision, and the next ridge is extracted with the energy recomputed.
Component \eqn{k} sums the original coefficients of its band,
\deqn{c_k[n] = \frac{1}{C} \mathrm{Re} \sum_{j = C_{k,n} - w}^{C_{k,n} + w} T_{n,j},}{c_k[n] = Re sum_{j = C_{k,n} - w}^{C_{k,n} + w} T_{n,j} / C,}
over the cells of the band, with
\eqn{C = C_\psi / (\sqrt{2\pi} \, 2 \ln 2)}{C = C_psi / (sqrt(2 pi) 2 ln 2)} and
\eqn{C_\psi = \int_0^\infty \hat\psi(\xi) / \xi \, d\xi}{C_psi = int_0^Inf psi^(xi) / xi dxi} for the Morlet
wavelet of center \code{morletCenter}. The bands may overlap, the
components need not add up to the signal, and no residue is returned. The
components are defined by the selected ridge bands, without imposing the
IMF conditions of EMD.
}
\value{
Double matrix, samples by components, one column per ridge.
}
\section{Conditions}{
A transform of zero energy stops the function with an error of class
\code{xemdZeroEnergy}. A sample where no cell satisfies the recurrence
within the tolerance stops it with an error.
}
\references{
Daubechies, I., Lu, J. and Wu, H.-T. (2009). Synchrosqueezed wavelet
transforms: a tool for empirical mode decomposition. arXiv:0912.2437.

Brevdo, E. and Thakur, G. \emph{Synchrosqueezing Toolbox}, revision
7e9fec0c6c9ed478dafac4479c0658d24fc6d8ef, functions curve_ext (original
code by Jianfeng Lu), curve_ext_multi and synsq_cwt_iw,
\url{https://github.com/ebrevdo/synchrosqueezing}.

The extraction of the ridges is adapted from curve_ext (original author
Jianfeng Lu, maintainer Eugene Brevdo) and curve_ext_multi of the
Synchrosqueezing Toolbox, as ported in SynchWave 1.1.2, revision
9ac80b2c3e60b417d0ee5567a66e8e64faa1a979 (port by Dongik Jang), and the
inversion from SynchWave 1.1.2 (MATLAB original by Eugene Brevdo; R port by
Dongik Jang, Hee-Seok Oh and Donghoh Kim), LGPL-2.0-or-later, which \code{xemd}
distributes under GPL-3 as section 3 of the GNU Library General Public
License version 2 permits; see the NOTICE file of the installed package.
}
\seealso{
\code{\link{transformSST}}, \code{\link{decomposeEWT}}
}
\examples{
SamplingInterval <- 1 / 1000
Time <- seq_len(1000) * SamplingInterval
Signal <- cos(2 * pi * 2 * Time) + 0.25 * cos(2 * pi * 24 * Time) +
    0.0625 * cos(2 * pi * 288 * Time)
SST <- transformSST(Signal, samplingInterval = SamplingInterval,
                    voiceCount = 16, morletCenter = 2 * pi,
                    threshold = sqrt(.Machine$double.eps))
Modes <- extractSSTModes(SST$transform, frequency = SST$frequency,
                         modeCount = 3, reconstructionWidth = 0,
                         morletCenter = 2 * pi)
dim(Modes)
sqrt(colMeans(Modes^2))
}
