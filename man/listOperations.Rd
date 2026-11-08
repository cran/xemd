\name{listOperations}
\alias{listOperations}
\title{Exported Operations and Their Arguments}
\description{
Lists the exported functions of the package with the names of their formal
arguments.
}
\usage{
listOperations()
}
\details{
The list reads the exports of the loaded namespace; it performs no
decomposition and reads no files.
}
\value{
A list with one element per exported function, in the order of
\code{sort()} on their names; each element is a list with \code{name}, the
function name, and \code{arguments}, the names of its formal arguments.
}
\seealso{
\code{\link{describeModule}}, \code{\link{listDecompositions}}
}
\examples{
Operations <- listOperations()
length(Operations)
Operations[[1]]
}
