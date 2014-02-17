\name{smartsSearchOB}
\alias{smartsSearchOB}
\title{
	SMARTS Search OB
}
\description{
	Perform searches for SMARTS patterns using Open Babel (requires
	ChemmineOB package to be installed).
}
\usage{
smartsSearchOB(sdfset, smartsPattern, uniqueMatches = TRUE)
}
\arguments{
  \item{sdfset}{
	  An SDFset of the compounds you want to search
}
  \item{smartsPattern}{
	  The SMARTS pattern as a string.
}
  \item{uniqueMatches}{
	  If true, only return the number of distinct matches, otherwise
	  return the number of all matches. 
}
}
\value{
	Returns a vector of counts, one for each input compound. 
}
\author{
	Kevin Horan
}


\examples{
	\dontrun{
		library(ChemmineOB)
		data(sdfsample)
		#look for rotable bonds
		rotableBonds = smartsSearchOB(sdfsample[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)
	}
}
