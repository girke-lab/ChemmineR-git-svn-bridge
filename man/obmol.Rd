\name{obmol}
\alias{obmol}
\alias{obmol-methods}
\title{
	obmol
}
\description{
Return reference to an OBMol from OpenBabel, if available. Operates
on SDF or SDFset objects.
}
\usage{
obmol(x)
}
\arguments{
  \item{x}{
object of class \code{SDF} or \code{SDFset}
}
}
\value{
	A pointer to an OBMol object, or a vector of pointers for an
	\code{SDFset}.
}
\author{
	Kevin Horan
}

\seealso{
\code{header}, \code{atomcount}, \code{bondblock}, \code{datablock}, \code{cid}, \code{sdfid} 
}
\examples{
## SDF/SDFset instances
if(require(ChemmineOB)){
	data(sdfsample)
	sdfset <- sdfsample
	sdf <- sdfset[[1]]

	obmolRef = obmol(sdf)
}
}

\keyword{ utilities }
