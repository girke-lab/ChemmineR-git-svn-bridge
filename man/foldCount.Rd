\name{foldCount}
\alias{foldCount}
\title{
	foldCount
}
\description{
	Returns the number of times this fingerprint has been folded.
}
\usage{
foldCount(x)
}
\arguments{
  \item{x}{
		Either an \code{FP} or an \code{FPset} object.
	}
}
\value{
	Returns the number of times this fingerprint has been folded.
}
\author{
	Kevin Horan
}
\examples{
	fp = new("FP",fp=c(1,0,1,1, 0,0,1,0))
	foldedFp=fold(fp)
	fc = foldCount(foldedFp) # == 1
}
