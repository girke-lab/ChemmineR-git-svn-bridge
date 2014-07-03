\name{result}
\alias{result}
\title{
Obtain the resulting output data from a ChemMine Tools Job
}
\description{
Accepts a \code{jobToken} job as returned by the function \code{launchCMTool} and returns
the final result. If the job is still running, the function will loop until the job is ready.
}
\usage{
result(object)
}
\arguments{
  \item{object}{
A \code{jobToken} job as returned by the function \code{launchCMTool}
}
}
\value{
Output will be in the format specified for this tool, as listed with the \code{listCMTools} function.
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{toolDetails}, \code{listCMTools}, \code{launchCMTool}, \code{browseJob}, \code{status}
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
