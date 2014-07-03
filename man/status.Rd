\name{status}
\alias{status}
\title{
Get Status of a ChemMine Tools Job
}
\description{
Returns the status of a launched ChemMine Tools job as represented by a \code{jobToken} object.
}
\usage{
status(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{jobToken} job as returned by the function \code{launchCMTool}
}
}
\value{
The status of the specified job is returned as a string. Possible values include "RUNNING", "FINISHED", or "FAILED".
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{toolDetails}, \code{listCMTools}, \code{launchCMTool}, \code{browseJob}, \code{result}
}
\examples{
\dontrun{
## list available tools
listCMTools()

## get detailed instructions on using a tool
toolDetails("Fingerprint Search")

## download compound 2244 from PubChem
job1 <- launchCMTool("pubchemID2SDF", 2244)

## check job status and download result
status(job1)
result1 <- result(job1)
}
}
\keyword{ utilities }
