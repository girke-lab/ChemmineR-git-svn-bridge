\name{propOB}
\alias{propOB}
\title{
	Properties from OpenBabel
}
\description{
	Generates the following descriptors in the given order:
	'abonds', 'atoms', 'bonds', 'dbonds', 'HBA1', 'HBA2', 'HBD', 'logP', 'MR', 'MW', 'nF', 'sbonds', 'tbonds', 'TPSA'.
}
\usage{
propOB(sdfSet)
}
\arguments{
  \item{sdfSet}{
	  An SDFset object.
	}
}
\value{
	A data frame with a row for each compound in the given data frame and a named
	column for each property.
}
\author{
	Kevin Horan
}

\examples{
	library(ChemmineR)
	data(sdfsample)
	propOB(sdfsample)
}
