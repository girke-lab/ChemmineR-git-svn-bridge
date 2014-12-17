\name{canonicalize}
\alias{canonicalize}
\title{
	Canonicalize
}
\description{
	Canonicalizes the atom numbering of a compound. The implimentation
	of this function is in Open Babel and requires the ChemmineOB
	package to function. 
}
\usage{
canonicalize(sdf)
}
\arguments{
  \item{sdf}{
	  Any sdfset object.
	}
}
\value{
	A new SDFset in which all compounds have been canonicalized
}
\references{
	http://openbabel.org/api/2.3/canonical_code_algorithm.shtml
}
\author{
	Kevin Horan
}
\seealso{
\code{\link{canonicalNumbering}}
}

\examples{
	data(sdfsample)
	canonicalSdf = canonicalize(sdfsample[1])
}
\keyword{canonical}
\keyword{morgan}
