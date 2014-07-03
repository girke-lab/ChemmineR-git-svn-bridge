\name{launchCMTool}
\alias{launchCMTool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Launch a Tool on ChemMine Tools
}
\description{
Accepts a tool name (string), input options, and input data to launch a remote
web tool on the ChemMine Tools website.
}
\usage{
launchCMTool(tool_name, input = "", ...)
}
\arguments{
  \item{tool_name}{
A tool name matching verbatim an existing tool name as listed by \code{listCMTools}.
}
  \item{input}{
Input data in the format required for this tool as listed by \code{listCMTools}.
}
  \item{\dots}{
Additional options as mentioned by running \code{toolDetails} on the tool specified.
}
}
\details{
By running the function \code{toolDetails} on a tool of choice, you can see
a pre-generated example function call for this tool.
}
\value{
\item{jobToken}{ for details see ?"jobToken-class"}
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{toolDetails}, \code{listCMTools}, \code{result}, \code{browseJob}, \code{status}
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
