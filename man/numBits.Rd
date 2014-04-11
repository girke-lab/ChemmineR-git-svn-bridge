\name{numBits}
\alias{numBits}
\title{
	numBits
}
\description{
	Returns the number of bits in a fingerprint.
}
\usage{
numBits(x)
}
\arguments{
  \item{x}{
		Either an \code{FP} or an \code{FPset} object.
	}

}
\value{
	The number of bits in this fingerprint object.
}
\author{
	Kevin Horan
}
\examples{
	fp = new("FP",fp=c(1,0,1,1, 0,0,1,0))
	n = numBits(fp) # == 8
}
