\name{decomposeTSL}
\alias{decomposeTSL}
\encoding{UTF-8}
\title{Decomposition of Time-Series Records}
\description{
Decomposes the series of a long table of records, grouped by record and
quantity, with one method of the package, and returns homogeneous tables for
all the methods.
}
\usage{
decomposeTSL(TSL, method, control, by)
}
\arguments{
\item{TSL}{Long table, a data frame or a \pkg{data.table}, with rows, unique
column names and the columns \code{t}, the time in seconds; \code{s}, the
sample; \code{ID}, the quantity, such as \code{"AT"}, \code{"VT"} or
\code{"DT"}; and \code{OCID}, the component. \code{t} and \code{s} are
finite numbers; \code{ID} and \code{OCID} are character without \code{NA}.
The other columns are metadata; they cannot take the names \code{mode},
\code{iterations}, \code{converged}, \code{stopReason},
\code{innerStopReason}, \code{frequency} or \code{boundary}.
Also reserved for the selected method are \code{realization} for EEMD;
\code{realization} and \code{sign} for CEEMD; \code{phase} for MaskedEMD;
and \code{direction} for BEMD. \code{TSL} is
not modified.}
\item{method}{One of \code{\link{listDecompositions}()}.}
\item{control}{Named list of arguments of \code{method}, common to all the
groups, except those the operation derives for each group:
\code{signal}, \code{samplingInterval}, \code{noise} and
\code{noiseChannels}. The ensembles take \code{realizationCount} in place
of \code{noise}; \code{\link{decomposeNAMEMD}} takes
\code{noiseChannelCount} and \code{noiseRatio}, a finite number
\eqn{\ge 0}, in place of \code{noiseChannels}. Controls in the units of the
signal are rejected, because AT, VT and DT do not share units: the rule
\code{"pair"} of \code{\link{decomposeMaskedEMD}} and any
\code{stopping$stopRule} of \code{\link{decomposeSVMD}} other than
\code{"modeCount"}.}
\item{by}{Character vector naming the metadata columns that identify a
record, or \code{character()} for a single record.}
}
\details{
A group is a record and an \code{ID}; its \code{OCID} are the channels of
the group and must share one time grid. Groups follow \code{by} ascending
and the IDs \code{"AT"}, \code{"VT"} and \code{"DT"} first, then the other
IDs ascending; the \code{OCID} of a group ascend. The order of the rows of
\code{TSL} does not intervene. A series needs two samples or more, a time
cannot repeat within a series, and the metadata other than \code{by} must be
constant within each record, \code{ID} and \code{OCID}.

The grid of a group is uniform when \code{diff(t)} equals the mean step
\eqn{\Delta t = (t_N - t_1)/(N - 1)} under \code{all.equal}; it is then used
unchanged, with \code{samplingInterval} \eqn{\Delta t}. A grid that is not
uniform is regularized with one warning per group: the step becomes
\eqn{1/F_s}, with \eqn{F_s} the smallest \eqn{2^a 5^b} not below
\eqn{1/\Delta t}; the grid is \eqn{t_1 + k/F_s}, and each \code{OCID} is
interpolated on it by \code{stats::splinefun(method = "monoH.FC")}. The
table \code{signal} holds the series actually decomposed.

Scalar methods decompose each \code{OCID} separately with the same
controls. The coupled methods decompose the matrix of the \code{OCID} of a
group in ascending order, without standardization:
\code{\link{decomposeBEMD}} needs exactly two \code{OCID} and
\code{\link{decomposeMEMD}} at least two; their \code{iterations},
\code{converged}, \code{stopReason} and \code{centers} are common and repeat
for every \code{OCID}. The ensembles draw \code{realizationCount} unit
normal realizations for each \code{OCID} from the random stream of the
caller, in the order group, \code{OCID}, realization.
\code{\link{decomposeNAMEMD}} draws \code{noiseChannelCount} channels for
each group and scales them by \code{noiseRatio} times
\eqn{\sqrt{\mathrm{mean}(s_c^2)}}{sqrt(mean(s_c^2))}, \eqn{s_c} the sample deviations of the
\code{OCID} of the group; the noise channels give no rows. A control valid
for one group and not for another, such as a center above the Nyquist
frequency of a group, stops the operation. Units are not inferred,
validated or converted; the metadata are carried to the output.

Channels are coupled within a record and an \code{ID}. The operation does
not impose a kinematic relationship between separate AT, VT and DT groups.
Equal mode indices across those groups do not establish corresponding
physical contributions.
}
\value{
A named list of \pkg{data.table} tables, each with the \code{by} columns,
the other metadata, \code{OCID} and \code{ID}, followed by
\item{signal}{\code{t}, \code{s}: the series decomposed, on the output
grid.}
\item{modes}{\code{mode}, an integer from 1, \code{t}, \code{s}. Each
\code{OCID} keeps its own number of modes, without padding; an \code{OCID}
without modes has no rows.}
\item{residue, fitResidual, meanNoise}{\code{t}, \code{s}, for the fields
that the method returns.}
\item{diagnostics}{\code{mode}, \code{iterations}, \code{converged}, when
the method returns iterations.}
\item{termination}{\code{stopReason}, and \code{innerStopReason} for
\code{\link{decomposeIteratedMaskedEMD}}, when the method returns a
termination. EEMD, CEEMD, MaskedEMD and BEMD carry the columns and rows of
their native \code{termination} data frame, with the series keys added;
the BEMD row repeats for both channels. In these rows \code{mode} identifies
the next attempted extraction, not a returned mode; \code{diagnostics}
continues to describe only returned modes.}
\item{centers, maskFrequencies}{\code{mode}, \code{frequency} in Hz.}
\item{boundaries}{\code{boundary}, \code{frequency} in Hz.}
Rows follow the \code{by} columns, \code{OCID}, \code{ID}, \code{mode} or
\code{boundary}, and \code{t}.
}
\references{
The regularization of the grid is adapted from \code{.regularize} and
\code{.rationalize} of gmsp 0.8.0, copyright (c) 2026
Alejandro Verri Kozlowski, MIT License: gmsp starts from the smallest
difference of the times and returns the grid from zero, whereas here the
step starts from the mean and the grid keeps the origin \eqn{t_1}. See the
NOTICE file of the installed package.
}
\section{Conditions}{
A group whose channels do not share one time grid, or whose series has
fewer than two samples, stops the function with an error of class
\code{xemdTimeGrid}. The warning about a non-uniform time grid is of class
\code{xemdWarning}. Both conditions name the group in their \code{unit}
field.
}
\seealso{
\code{\link{listDecompositions}}, \code{\link{decomposeEMD}},
\code{\link{decomposeMVMD}}
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
OUT <- decomposeTSL(TSL, method = "decomposeEMD",
                    control = list(sifting = list(stopRule = "doubleThreshold")),
                    by = character())
names(OUT)
OUT$termination
head(OUT$modes)

# A coupled method: three channels of one record decomposed together.
SamplingInterval <- 1 / 1000
TimeV <- seq_len(1000) * SamplingInterval
Series <- cbind(cos(2 * pi * 2 * TimeV) + 0.25 * cos(2 * pi * 24 * TimeV) +
                    0.0625 * cos(2 * pi * 288 * TimeV),
                cos(2 * pi * 2 * TimeV), 0.25 * cos(2 * pi * 24 * TimeV))
Record <- data.table::data.table(ID = "AT",
                                 OCID = rep(c("H1", "H2", "UP"), each = 1000),
                                 t = rep(TimeV, 3), s = as.vector(Series))
Coupled <- decomposeTSL(Record, method = "decomposeMVMD",
                        control = list(initialCenters = (0:2) / (6 * SamplingInterval),
                                       alpha = 2000 * SamplingInterval^2 / (8 * pi^2),
                                       tau = 0.5, fixedZeroCenter = FALSE,
                                       tolerance = 1e-6, maxIterations = 98),
                        by = character())
Coupled$centers
}
