\name{measureModesTSL}
\alias{measureModesTSL}
\encoding{UTF-8}
\title{Metrics of the Modes of Time-Series Records}
\description{
Measures every mode and every series of the tables of
\code{\link{decomposeTSL}}, for any of its methods: frequency, bandwidth,
energy, energy fraction and mean period of each mode, and number of modes
and energy fractions of the remainder and the energy difference described
by \code{crossFraction} for each series.
}
\usage{
measureModesTSL(decomposition)
}
\arguments{
\item{decomposition}{The list of tables of \code{\link{decomposeTSL}}, with
\code{signal}, \code{modes} and at most one of \code{residue} and
\code{fitResidual}. Its keys are the columns of \code{signal} other than
\code{t} and \code{s}: the \code{by} columns, the other metadata,
\code{OCID} and \code{ID}; they cannot take the names \code{mode},
\code{frequency}, \code{bandwidth}, \code{energy}, \code{energyFraction},
\code{meanPeriod}, \code{modeCount}, \code{residueFraction} or
\code{crossFraction}. \code{t} and \code{s} are finite numbers and
\code{mode} holds whole numbers from 1. Each series of \code{signal} has two
samples or more at distinct times, and each mode and remainder of a series
holds exactly the times of that series. \code{decomposition} is not
modified.}
}
\details{
A mode \eqn{c} of a series \eqn{x} with \eqn{N} samples uses the step
\eqn{\Delta t = (t_N - t_1)/(N - 1)} of the grid of its group, with its
samples ordered by \code{t}; the order of the rows does not intervene. Its
spectrum is that of the VMD family: the mode is extended to \eqn{2N}
samples by reflecting \eqn{\lceil N/2 \rceil}{ceiling(N/2)} samples at the
start and the others at the end, \eqn{P_j} is the squared modulus of the
unnormalized discrete Fourier transform of the extension and
\eqn{|\omega_j| = 2\pi \min(j, 2N - j)/(2N \Delta t)}{|omega_j| = 2 pi
min(j, 2N - j) / (2N dt)}, \eqn{j = 0, \ldots, 2N - 1}. Then
\code{frequency} is the power centroid
\deqn{f = \frac{\sum_j |\omega_j| P_j}{2\pi \sum_j P_j}}{f = sum_j
|omega_j| P_j / (2 pi sum_j P_j)}
and \code{bandwidth} the power deviation around it on the same grid,
\deqn{b = \left[\frac{\sum_j (|\omega_j|/2\pi - f)^2 P_j}{\sum_j
P_j}\right]^{1/2},}{b = [sum_j (|omega_j| / (2 pi) - f)^2 P_j / sum_j
P_j]^(1/2),}
both in Hz. \code{frequency} equals the \code{centers} of
\code{\link{decomposeVMD}} and \code{\link{decomposeSVMD}}, except a center
fixed at zero, which is not a centroid; in \code{\link{decomposeMVMD}} it is
the centroid of each channel, not the common center. A mode of zero energy
determines no centroid: its \code{frequency} and \code{bandwidth} are
\code{NA}. \code{energy} is \eqn{\Delta t \sum_n c[n]^2}{dt sum_n c[n]^2},
\code{energyFraction} is \eqn{\sum c^2 / \sum x^2}{sum(c^2) / sum(x^2)}, and
\code{meanPeriod} is \eqn{N \Delta t}{N dt} over the number of maxima of
the mode, or \eqn{N \Delta t}{N dt} without maxima.

For each series, \code{modeCount} is the number of its modes;
\code{residueFraction} is \eqn{\sum r^2 / \sum x^2}{sum(r^2) / sum(x^2)},
with \eqn{r} its \code{residue} or \code{fitResidual}, and zero when the
method returns neither, as \code{\link{decomposeEWT}}; and
\code{crossFraction} is
\deqn{1 - \frac{\sum_k \sum c_k^2 + \sum r^2}{\sum x^2},}{1 - (sum_k
sum(c_k^2) + sum(r^2)) / sum(x^2),}
the input-energy fraction left after subtracting modal and remainder
energies. It equals the normalized cross products when the modes and
remainder reconstruct that input. Otherwise it also reflects the
reconstruction discrepancy, as with EEMD's retained \code{meanNoise}.
A series without modes has \code{modeCount} 0, and its remainder is the
series.

Energy here is a quadratic signal measure, with squared signal units times
seconds. Its identification with physical energy requires a relation for
the measured quantity. A large modal energy fraction does not establish
that the component is free of contamination.
}
\value{
A named list of two \pkg{data.table} tables, each with the keys of
\code{decomposition} followed by
\item{metrics}{\code{mode}, \code{frequency} in Hz, \code{bandwidth} in Hz,
\code{energy} in the squared units of the series times seconds,
\code{energyFraction} and \code{meanPeriod} in seconds; one row per mode,
in the order of \code{modes}.}
\item{metricsSeries}{\code{modeCount}, \code{residueFraction} and
\code{crossFraction}; one row per series, in the order of \code{signal}.}
}
\section{Conditions}{
A series of \code{signal} with zero energy leaves its energy fractions
undefined and stops the function with an error of class
\code{xemdZeroEnergy} whose \code{unit} field names the series.
}
\seealso{
\code{\link{decomposeTSL}}, \code{\link{selectModesTSL}},
\code{\link{reconstructTSL}}
}
\examples{
Time <- seq(0, 9, length.out = 3000)
Signal <- sin(pi * Time) + sin(2 * pi * Time) + sin(6 * pi * Time) +
    0.5 * Time
# The same series under three IDs and three components.
TSL <- data.table::rbindlist(lapply(c("AT", "VT", "DT"), function(id) {
    data.table::rbindlist(lapply(c("H1", "H2", "UP"), function(ocid) {
        data.table::data.table(ID = id, OCID = ocid, t = Time, s = Signal)
    }))
}))
Sifting <- list(stopRule = "doubleThreshold")
Decomposition <- decomposeTSL(TSL, method = "decomposeEMD",
                              control = list(sifting = Sifting),
                              by = character())
OUT <- measureModesTSL(Decomposition)
OUT$metrics[OUT$metrics$OCID == "H1" & OUT$metrics$ID == "AT", ]
OUT$metricsSeries
}
