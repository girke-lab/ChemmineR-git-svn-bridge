\name{jarvisPatrick_c}
\alias{jarvisPatrick_c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Jarvis Patrick Clustering in C code
}
\description{
	This not meant to be used directly, use \code{\link{jarvisPatrick}} instead. It is exposed so other
	libraries can make use of it.
}
\usage{
jarvisPatrick_c(neighbors,minNbrs,fast=TRUE,bothDirections=FALSE,linkage = "single")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{neighbors}{
	  A matrix of integers. Non integer matricies will be coerced. 
	  Each row represensts one element, indexed 1 to N. The values in row
	  i should be the index value of the neighbors of i. Thus, each value should 
	  itself be a valid row index. 
}
  \item{minNbrs}{
	  The minimum number of common neibhbors needed for two elements to be merged.
}
  \item{fast}{
	  If true, only the neibhors given in each row are checked to see if they share \code{minNbrs} 
	  neighbors in common. If false, all pairs of elements are compared. For a matrix of size NxM,
	  the first method yeilds a running time of O(NM), while the second yeilds a running time of O(N^2).
}
  \item{bothDirections}{
	  If true, two elements must contain each other in their neighbor list in order to be merged.
	  If false and fast is true, then only one element must contain the other as a neighbor. If
	  false and fast is false, than neither element must contain the other as a neighbor, though
	  in all cases there must still be at least \code{minNbrs} neibhros in common.
}
	\item{linkage}{
		See \code{\link{jarvisPatrick}} for details.
	}
}
\value{
	A cluster array with no names.
}
\author{
	Kevin Horan
}
