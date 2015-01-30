\name{setPriorities}
\alias{setPriorities}
\alias{forestSizePriorities}
\alias{randomPriorities}
\title{
	Set Priorities
}
\description{
	This function should be run after loading a complete set of data.
	It will find each group of compounds which share the same
	descriptor and call the given function, \code{priorityFn}, 
	with the compound_id numbers of the group. This function should
	then assign priorities to each compound-descriptor pair, however
	it wishes. Priorities are integer values with lower values being
	used in preference of higher values.

	It is important that this function be called after all data is
	loaded. It may be that a compound loaded at the beginning of a data
	set shares a descriptor with a compound loaded near the end of the
	data set. If the priorities were set at some point in between these
	then it would not see all the compounds for that one descriptor.

	If a SNOW cluster and connection source function are given, it will
	run in parallel.

	Some pre-defined functions that can be use for \code{priorityFn}
	are:

	\code{randomPriorities}: Set the priorities of compounds within a descriptor group
	randomly.

	\code{forestSizePriorities}: Set the priority based on the number
	of disconnected components (trees) within the compound. Compounds
	with fewer trees will have a higher priority (lower numerical
	value) than compounds with more trees.
}

\usage{
setPriorities(conn,priorityFn,descriptorIds=c(),cl=NULL,connSource=NULL){
forestSizePriorities(conn,compIds)
randomPriorities(conn,compIds)
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
	\item{descriptorIds}{
			If given then only re-compute priorities for groups
			involving descriptors in this list. This is useful
			for updating priorities after adding new compounds
			to an existing database.
	}
	\item{cl}{
		A SNOW cluster on which to run jobs on.
	}
	\item{connSource}{
		A function to create a new database connection with. This will be run once for each new job
		created. It must return a newly created connection, not a reference to an existing connection.
	}
	\item{compIds}{
	  The compound_id values for each group.
	}

}
\value{
	For \code{setPriorities}, no value is returned.
	\code{randomPriorities} and \code{forestSizePriorities} return
	a data.frame with columns "compound_id" and "priority".
}
\author{
	Kevin Horan
}
\examples{
	\dontrun{
		data(sdfsample)
		conn = initDb("sample.db")
		sdfLoad(conn,sdfsample)
		setPriorities(conn,forestSizePriorities)
	}
}
