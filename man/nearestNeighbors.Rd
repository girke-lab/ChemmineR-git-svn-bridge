\name{nearestNeighbors}
\alias{nearestNeighbors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Nearest Neighbors
}
\description{
   Computes the nearest neighbors of descriptors in an FPset or APset object for use with the \code{\link{jarvisPatrick}} clustering
   function. Only one of \code{numNbrs} or \code{cutoff} should be given, \code{cutoff} will take precedence if
   both are given. If \code{numNbrs} is given, then that many neighbors will be returned for each item in the set.
   If \code{cutoff} is given, then, for each item X, every neighbor that has a similarity value greater than or equal to 
   the cutoff will be returned in the neighbor list for X.
}
\usage{
   nearestNeighbors(x, numNbrs = NULL, cutoff = NULL, ...)
}
\arguments{
  \item{x}{
     Either an FPset or an APset.
   }
  \item{numNbrs}{
     Number of neighbors to find for each item. If not enough neighbors can be found the 
     matrix will be padded with NA.
   }
  \item{cutoff}{
      The minimum similarity value an item must have to another item in order to be included in that
      items neighbor list. This parameter takes precedence over \code{numNbrs}.
      This parameter allows to obtain tighter clustering results.
   }
  \item{\dots}{
     These parameters will be passed into the distance function used, either \code{cmp.similarity} or 
     \code{fpSim}, for APset and FPset, respectively.
   }
}
\value{
   The return value is a list with the following components:
   \item{ids}{index values of nearest neighbors, for each item. If \code{cutoff} is
      used, this will be a list of lists, otherwise it will be a matrix}
   \item{names}{The names of each item in the set, as returned by cid}
   \item{similarities}{The similarity values of each neighbor to the item for that row. 
      This will also be either a list of lists or a matrix, depending on whether or not
      \code{cutoff} was used. Each similarity values corresponds to the id number in the 
      same position in the ids entry}
}
\author{
   Kevin Horan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
   \code{\link{jarvisPatrick}}
   \code{\link{trimNeighbors}}
}
\examples{
   data(sdfsample)
   ap = sdf2ap(sdfsample)
   nnm = nearestNeighbors(ap,cutoff=0.5)
   clustering = jarvisPatrick(nnm,k=2,mode="a1b")
}
