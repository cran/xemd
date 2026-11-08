\name{listDecompositions}
\alias{listDecompositions}
\title{Methods of the Record Operation}
\description{
Returns the names of the decomposition methods that \code{\link{decomposeTSL}}
accepts.
}
\usage{
listDecompositions()
}
\details{
The scalar decompositions are applied to each \code{OCID} of a group and
the coupled ones (\code{\link{decomposeBEMD}}, \code{\link{decomposeMEMD}},
\code{\link{decomposeNAMEMD}} and \code{\link{decomposeMVMD}}) to the
\code{OCID} of a group together. \code{\link{transformSST}}, the analysis
operations, \code{\link{resiftModes}} and
\code{\link{estimateMaskFrequency}} have no record operation.
}
\value{
Character vector of function names: the thirteen scalar methods followed by
the four coupled ones.
}
\seealso{
\code{\link{decomposeTSL}}, \code{\link{listOperations}}
}
\examples{
listDecompositions()
}
