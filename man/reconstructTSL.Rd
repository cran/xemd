\name{reconstructTSL}
\alias{reconstructTSL}
\encoding{UTF-8}
\title{Reconstruction of Time-Series Records from Their Modes}
\description{
Adds up, for each series of the tables of \code{\link{decomposeTSL}}, the
modes that a selection retains and, unless excluded, the remainder of the
series.
}
\usage{
reconstructTSL(decomposition, selection, includeResidue = TRUE)
}
\arguments{
\item{decomposition}{The list of tables of \code{\link{decomposeTSL}}, as
for \code{\link{measureModesTSL}}. It is not modified.}
\item{selection}{A table with the keys and \code{mode} of every mode of
\code{decomposition}, once each, and \code{retained}, \code{TRUE} or
\code{FALSE}, as the table \code{selection} of
\code{\link{selectModesTSL}}. It is not modified.}
\item{includeResidue}{\code{TRUE}, the default, adds the remainder of each
series, \code{residue} or \code{fitResidual}; \code{FALSE} leaves it out.}
}
\details{
The reconstruction of a series is the sum of its retained modes and, with
\code{includeResidue}, its remainder. A series without retained modes sums
zero, and a method without remainder, as \code{\link{decomposeEWT}}, adds
none. With every mode and the remainder, the result follows the
reconstruction relation of the selected method. For EMD the sum recovers
the input up to rounding. In \code{\link{decomposeEEMD}} with the divisor
\code{"realizations"}, it gives the signal plus \code{meanNoise}.
The order of the rows of the tables does not intervene.

A partial sum implements the declared selection; neither the sum nor the
choice to include a remainder establishes that the result is free of noise.
}
\value{
A named list with one \pkg{data.table} table, \code{reconstruction}: the
keys of \code{decomposition}, \code{t} and \code{s}, on the grid and in the
order of \code{signal}.
}
\seealso{
\code{\link{selectModesTSL}}, \code{\link{measureModesTSL}},
\code{\link{decomposeTSL}}
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
Metrics <- measureModesTSL(Decomposition)$metrics
# A minimum energy fraction of zero retains every mode: with the residue,
# the reconstruction is the signal.
Every <- selectModesTSL(Decomposition, metrics = Metrics,
                        criterion = list(rule = "energyFraction", minimum = 0,
                                         action = "keep"))
Complete <- reconstructTSL(Decomposition, selection = Every$selection)
max(abs(Complete$reconstruction$s - Decomposition$signal$s))
# The modes in [2, 15] Hz, without the residue.
Band <- selectModesTSL(Decomposition, metrics = Metrics,
                       criterion = list(rule = "band", fMin = 2, fMax = 15,
                                        action = "keep"))
Filtered <- reconstructTSL(Decomposition, selection = Band$selection,
                           includeResidue = FALSE)
head(Filtered$reconstruction)
}
