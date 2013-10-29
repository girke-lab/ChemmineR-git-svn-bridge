\name{maximallyDissimilar}
\alias{maximallyDissimilar}
\title{
	Maximally Dissimilar
}
\description{
	Find a set of compounds that are far away from each other.
}
\usage{
maximallyDissimilar(compounds, n, similarity = cmp.similarity)
}
\arguments{
  \item{compounds}{
	  The set of items from which to pick \code{n} dissimlar items. This can be a list of
	  anything that the similarity function will accept. By default this will be
	  an APset.
	}
  \item{n}{
	  The number of dissimilar items to return.
	}
  \item{similarity}{
	  The similarity function to use. By default Tanimoto will be used on APset objects.
	  Internally, this will be converted to a distance function using \code{1-similarity(a,b)}, so
	  whatever similarity function you use should return a value between 0 and 1.
	}
}
\details{
	This will run in O(length(compounds)n) time. Based on the algorithm described in (Higgs,1997).
}
\value{
	A vector of indexes of the dissimilar items.
}
\references{
	Higgs, R.E., Bemis, K.G., Watson, I.A., and Wikel, J.H. 1997. Experimental designs for
	selecting molecules from large
	chemical databases. J. Chem. Inf. Comput. Sci. 37, 861-870

}
\author{
	Kevin Horan
}

\examples{
	data(apset)
	maximallyDissimilar(apset,10)
}
\keyword{dissimilar}
