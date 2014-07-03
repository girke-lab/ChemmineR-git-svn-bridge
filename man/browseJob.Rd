\name{browseJob}
\alias{browseJob}
\title{
Open ChemMine Tools Job in Web Browser
}
\description{
Launches a web browser to view the results of a ChemMine Tools web job with an interactive online viewer.
Note that this reassigns the job to the current logged in user within the browser, so it becomes no longer
accessible by the \code{result} and \code{status} functions. Any results should be saved within R before launching a browser.
}
\usage{
browseJob(object)
}
\arguments{
  \item{object}{
A \code{jobToken} job as returned by the function \code{launchCMTool}
}
}
\value{
Returns an URL string which can be used to access the job results. The function also attempts to open the url with the \code{browseURL} function. As this URL can only be used once, the returned string is only useful if the \code{browseURL} function fails to open a browser.
}
\references{
See ChemMine Tools at \url{http://chemmine.ucr.edu}.
}
\author{
Tyler William H Backman
}

\seealso{
Functions: \code{toolDetails}, \code{listCMTools}, \code{launchCMTool}, \code{result}, \code{status}
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

## open job in web browser
browseJob(job1)
}
}
\keyword{ utilities }
