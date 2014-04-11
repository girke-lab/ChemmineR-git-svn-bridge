\name{genAPDescriptors}
\alias{genAPDescriptors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Generate AP Descriptors
}
\description{
	Generates Atom Pair descriptors using a fast C function.
}
\usage{
genAPDescriptors(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
		A single SDF object.
	}
}
\value{
	A vector of descriptors for the compound given.  An AP object
	can be generated as shown in the example below.
}
\author{
	Kevin Horan
}
\examples{
	library(ChemmineR)
	data(sdfsample)
	sdf = sdfsample[[2]]
	ap = new("AP", AP=genAPDescriptors(sdf))
}
