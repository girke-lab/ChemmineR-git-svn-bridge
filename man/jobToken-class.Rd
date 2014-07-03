\name{jobToken-class}
\Rdversion{1.1}
\docType{class}
\alias{jobToken-class}
\alias{show,jobToken-method}

\title{Class \code{"jobToken"}}
\description{
Container for storing a reference to a remote job ran on the ChemMine Tools web server.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("jobToken", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{tool_name}:}{Object of class \code{"character"} }
    \item{\code{jobId}:}{Object of class \code{"character"} }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "jobToken")}: check the status of a launched job }
	 }
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{launchCMTool}, \code{toolDetails}, \code{listCMTools}, \code{result}, \code{browseJob}, \code{status}
}
\examples{
showClass("jobToken")
\dontrun{
## launch a job on the server and obtain jobToken back
job1 <- launchCMTool("pubchemID2SDF", 2244)

## check status of the job
status(job1)

## obtain results
result1 <- result(job1)
result1
}
}
\keyword{classes}
