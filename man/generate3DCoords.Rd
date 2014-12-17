\name{generate3DCoords}
\alias{generate3DCoords}
\title{
	Generate 3D Coords
}
\description{
	Uses Open Babel to compute 3D coordinates given an SDFset with
	only 2D coordinates.
}
\usage{
generate3DCoords(sdf)
}
\arguments{
  \item{sdf}{
	  Any sdfset object.
	}
}
\value{
	A new SDFset in which all compounds have 3D coordinates.
}
\author{
	Kevin Horan
}

\examples{
	data(sdfsample)
	sdf3D = generate3DCoords(sdfsample[1])
}
