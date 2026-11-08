\name{selectModesTSL}
\alias{selectModesTSL}
\encoding{UTF-8}
\title{Selection of the Modes of Time-Series Records}
\description{
Selects the modes of the tables of \code{\link{decomposeTSL}} with a
declared criterion: explicit mode indices, a minimum energy fraction or a
frequency band, to keep the modes that meet it or to remove them.
}
\usage{
selectModesTSL(decomposition, metrics, criterion)
}
\arguments{
\item{decomposition}{The list of tables of \code{\link{decomposeTSL}}, as
for \code{\link{measureModesTSL}}; its metadata cannot take the names
\code{mode}, \code{criterion} or \code{retained}. It is not modified.}
\item{metrics}{The table \code{metrics} of \code{\link{measureModesTSL}}
for the same decomposition: one row for each of its modes, each with an
\code{energyFraction}. It is not modified.}
\item{criterion}{A list with \code{rule}, exactly the fields of that rule,
and \code{action}, \code{"keep"} or \code{"remove"}. The rules are
\describe{
\item{\code{"modes"}}{with \code{modes}, a table with \code{mode}, whole
numbers from 1, and key columns of \code{decomposition}. Each row names the
mode of that index of every series whose keys it matches, so the indices go
by group (the \code{by} columns and \code{ID}) or by series (with
\code{OCID}); each row names a mode of the decomposition.}
\item{\code{"energyFraction"}}{with \code{minimum}, a finite number not
below zero: a mode meets it when its \code{energyFraction} is not below
\code{minimum}.}
\item{\code{"band"}}{with \code{fMin} and \code{fMax} in Hz,
\eqn{0 \le f_{\min} \le f_{\max}}{0 <= fMin <= fMax}: a mode meets the
closed band when \eqn{f_{\min} \le}{fMin <=} \code{frequency}
\eqn{\le f_{\max}}{<= fMax}. A mode without frequency, of zero energy, is in
no band.}
}
\code{"keep"} retains the modes that meet the criterion and
\code{"remove"} the others.}
}
\details{
The criterion is supplied by the caller. Mode indices, energy fractions and
centroid frequencies do not by themselves identify uncontaminated physical
components. Equal mode indices in independently decomposed series need not
represent the same physical contribution.

The band compares \code{frequency} alone, whereas the rules \code{keep.Fo}
and \code{remove.Fo} of gmsp 0.8.0 require the interval
\eqn{[F_o - BW, F_o + BW]}{[Fo - BW, Fo + BW]} inside the band.
}
\value{
A named list of two \pkg{data.table} tables:
\item{modes}{The rows of \code{modes} of \code{decomposition} of the
retained modes, in its order.}
\item{selection}{Every mode of \code{decomposition} with its keys,
\code{mode}, \code{criterion} and \code{retained}, in the order of the
decomposition. \code{criterion} is the criterion as declared text:
\code{"<action> listed modes"}, \code{"<action> energyFraction >= <minimum>"}
or \code{"<action> frequency in [<fMin>, <fMax>] Hz"}, such as
\code{"keep frequency in [2, 15] Hz"}.}
}
\seealso{
\code{\link{measureModesTSL}}, \code{\link{reconstructTSL}},
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
# Keep the modes whose frequency lies in [2, 15] Hz.
Band <- selectModesTSL(Decomposition, metrics = Metrics,
                       criterion = list(rule = "band", fMin = 2, fMax = 15,
                                        action = "keep"))
Band$selection[Band$selection$OCID == "H1" & Band$selection$ID == "AT", ]
# Remove the first mode of every series of the ID VT.
Listed <- data.frame(ID = "VT", mode = 1)
First <- selectModesTSL(Decomposition, metrics = Metrics,
                        criterion = list(rule = "modes", modes = Listed,
                                         action = "remove"))
table(First$selection$ID, First$selection$retained)
}
