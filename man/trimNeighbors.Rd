\name{trimNeighbors}
\alias{trimNeighbors}
\title{
   Trim Neighbors
}
\description{
   Further reduce the cutoff value of a nearest neighbor (NN) table, as produced by
   \code{\link{nearestNeighbors}}. This allows one to compute a very relaxed NN table
   initially, and then quickly restrict it later without having to re-compute all the
   similarities.
}
\usage{
   trimNeighbors(nnm, cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nnm}{
      A nearest neighbor table, as produced by \code{\link{nearestNeighbors}}.
   }
  \item{cutoff}{
      The new similarities cutoff value. All pairs with a similarity less than this value
      will be removed from the table.
   }
}
\value{
   The return value has the same structure as \code{nnm}, with some neighbors
   removed from the \code{indexes} and \code{similarties} entries.
}
\author{
   Kevin Horan
}

\seealso{
   \code{\link{jarvisPatrick}}
   \code{\link{nearestNeighbors}}
}
\examples{

   data(sdfsample)
   ap = sdf2ap(sdfsample)
   nnm = nearestNeighbors(ap,numNbrs=20)
   nnm = trimNeighbors(nnm,cutoff=0.5)
   clustering = jarvisPatrick(nnm,k=2,mode="a1b")
}
