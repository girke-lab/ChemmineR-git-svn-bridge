\name{randomPriorities}
\alias{randomPriorities}
\title{
	Random Priorities
}
\description{
	Set the priorities of compounds within a descriptor group
	randomly.

	This function is meant to be passed into
	\code{\link{setPriorities}}, it is not meant to be called directly.
}
\usage{
randomPriorities(compIds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compIds}{
	  The compound_id values for each group.
}
}
\value{
	A data.frame with columns "compound_id" and "priority".
}
\author{
	Kevin Horan
}

\seealso{
	\code{\link{setPriorities}}
}
\examples{
	\dontrun{
		data(sdfsample)
		conn = initDb("sample.db")
		sdfLoad(conn,sdfsample)
		setPriorities(conn,randomPriorities)
	}
}
