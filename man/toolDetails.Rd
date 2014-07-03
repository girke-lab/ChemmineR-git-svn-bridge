\name{toolDetails}
\alias{toolDetails}
\title{
Detailed instructions for each ChemMine Tools web tool
}
\description{
Connects to the ChemMine Tools web service, and provides detailed instructions and example
function calls for any tool.
}
\usage{
toolDetails(tool_name)
}
\arguments{
  \item{tool_name}{
A tool name matching verbatim an existing tool name as listed by \code{listCMTools}.
}
}
\details{
Prints instructions to console.
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{launchCMTool}, \code{listCMTools}, \code{result}, \code{browseJob}, \code{status}
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
