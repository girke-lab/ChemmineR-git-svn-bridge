\name{byCluster}
\alias{byCluster}
\title{
	By Cluster
}
\description{
	Re-organize a vector valued clustering into an list which groups cluster members together
}
\usage{
byCluster(clustering, excludeSingletons = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
		A named vector in which the names are cluster members and the values are cluster labels. This
		is format output by jarvisPatrick.
	}
  \item{excludeSingletons}{
	  If true only clusters with more than 1 member will be in the output, otherwise all clusters will be used.
	}
}
\value{
	A list with a slot for each cluster. Each slot of the list is a vector containing the cluster members.
}
\author{
	Kevin Horan
}


\seealso{
	\code{\link{jarvisPatrick}}
}
\examples{
	data(apset)
	cl = jarvisPatrick(nearestNeighbors(apset,cutoff=0.6),k=2)
	print(byCluster(cl))
}
