\name{exactMassOB}
\alias{exactMassOB}
\title{
	Exact Mass (Monoisotopic Mass)
}
\description{
	Computes the exact mass of each compound given.
}
\usage{
exactMassOB(sdfset)
}
\arguments{
  \item{sdfset}{
		Any SDFset object.
	}
}
\value{
	A vector of mass values.
}
\author{
	Kevin Horan
}
\examples{
	library(ChemmineR)
	data(sdfsample)
	mass = exactMassOB(sdfsample)
}
\keyword{exact}
\keyword{mass}
\keyword{monoisotopic}
