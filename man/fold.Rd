\name{fold}
\alias{fold}
\title{
	Fold
}
\description{
	Fold a fingerprint. This takes the second half of the fingerprints
	and combines with the first half with a logical 'OR' operation. The
	result is a fingerprint with half as many bits.
}
\usage{
fold(x, count = 1, bits = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The fingerprint(s) to fold. This can be either an \code{FP}
	  or an \code{FPset} object.
	}
  \item{count}{
	  The number of times to fold this fingerprint. Folding will stop
	  early if the fingerprint is reduced down to 1 bit before reaching
	  the requested fold count.
	}
  \item{bits}{
	  Fold this fingerprint until it is \code{bits} bits long. An
	  exception will be thrown if \code{bits} is not reachable.
	}
}
\value{
	The new, folded, fingerprint.
}
\author{
	Kevin Horan
}
\examples{
	fp = new("FP",fp=c(1,0,1,1, 0,0,1,0))
	foldedFp = fold(fp,bits=4)
}
