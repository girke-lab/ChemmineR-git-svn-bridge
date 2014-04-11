\name{fromNNMatrix}
\alias{fromNNMatrix}
\title{
	From Nearest Neighbor Matrix
}
\description{
	Converts a nearest neighbor matrix into a list that can be used with the jarvisPatrick function.
}
\usage{
fromNNMatrix(data, names = rownames(data))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
		A matrix containing integer valued indexes which represent items to be clustered. The index values contained in
		the matrix must be smaller than the number of rows in the matrix.
		Each row in the matrix represents one item and the columns are the nearest neighbors of that item.
	}
  \item{names}{
	  The names for each row. The rownames of data will be used if not given.
	}
}
\value{
	A list containing the slots "indexes" and "names".
}
\author{
	Kevin Horan
}

\seealso{
	\code{\link{jarvisPatrick}}
}
\examples{


	data(apset)

	nn = nearestNeighbors(apset,cutoff=0.6)
	nnMatrix = nn$indexes

	cl = jarvisPatrick(fromNNMatrix(nnMatrix),k=2)
}
