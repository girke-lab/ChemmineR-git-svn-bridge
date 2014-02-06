\name{setPriorities}
\alias{setPriorities}
\title{
	Set Priorities
}
\description{
	This function should be run after loading a complete set of data.
	It will find each group of compounds which share the same
	descriptor and call the given function, \code{priorityFn}, 
	with the compound_id numbers of the group. This function should
	then assign priorities to each compound-descriptor pair, however
	it wishes.

	It is important that this function be called after all data is
	loaded. It may be that a compound loaded at the beginning of a data
	set shares a descriptor with a compound loaded near the end of the
	data set. If the priorities were set at some point in between these
	then it would not see all the compounds for that one descriptor.
}

\usage{
setPriorities(conn, priorityFn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	  A database connection object.
}
  \item{priorityFn}{
	  This function will be called with the compound_id numbers 
	  associated with the same descriptor. It should use the id
	  numbers to lookup whatever data it wants to assign a priority to
	  each compound. These priority values will be used to pick a
	  compound to represent the group in cases where only one compound
	  is needed for each descriptor. 

	  The function should return a data.frame with the fields
	  "compound_id" and "priority". The order of the rows is not
	  important.
}
}
\value{
	No value is returned
}
\author{
	Kevin Horan
}
\seealso{
	\code{\link{randomPriorities}}
}
\examples{
	\dontrun{
		data(sdfsample)
		conn = initDb("sample.db")
		sdfLoad(conn,sdfsample)
		setPriorities(conn,randomPriorities)
	}
}
