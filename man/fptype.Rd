\name{fptype}
\alias{fptype}
\title{
	fptype
}
\description{
	Returns the type label of this fingerprint
}
\usage{
fptype(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Either an \code{FP} or an \code{FPset} object.
	}
}
\value{
	The type label of this fingerprint.
}
\author{
	Kevin Horan
}
\examples{

	fp = new("FP",fp=c(1,0,1,1, 0,0,1,0),type="testFP")
	type = fptype(fp) # == "testFP"
}
