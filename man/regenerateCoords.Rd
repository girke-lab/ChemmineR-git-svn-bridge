\name{regenerateCoords}
\alias{regenerateCoords}
\title{
	Re-generate 2D Coordinates
}
\description{
	This uses Open Babel (requires ChemmineOB package) to re-generate
	the 2D coordinates of compounds. This often results in a 
	nicer layout of the compound when plotting.
}
\usage{
regenerateCoords(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
	  A SDF or SDFset object whose coordinates will be re-generated.
}
}
\value{
	Either an SDF object if given an SDF, or else and SDFset.
}
\author{
	Kevin Horan
}

\seealso{
	\code{\link{plotStruc}}
}
\examples{
	\dontrun{
		data(sdfsample)
		prettySdfset = regenerateCoords(sdfsample[1:4])
	}

}
