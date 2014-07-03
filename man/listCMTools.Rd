\name{listCMTools}
\alias{listCMTools}
\title{
List all available ChemMine Tools
}
\description{
Connects to the ChemMine Tools web service and obtains a list of all available tools, and their
input and output formats.
}
\usage{
listCMTools()
}
\value{
\item{data.frame}{ A four column data.frame which describes a tool on each row}
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{toolDetails}, \code{launchCMTool}, \code{result}, \code{browseJob}, \code{status}
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
