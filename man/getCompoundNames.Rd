\name{getCompoundNames}
\alias{getCompoundNames}
\title{
	Get Compound Names
}
\description{
	Fetch the names of the given compound ids, if they exist
}
\usage{
getCompoundNames(conn, compoundIds, keepOrder = FALSE, allowMissing = FALSE)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{compoundIds}{
     A vector of compound ids.
   }
  \item{keepOrder}{
     If true, the order of the output compound ids will be the same as the 
     input names. This imposes a performance hit that can be significant for
     large datasets, thus it should be left FALSE unless needed.
   }
	\item{allowMissing}{
		When this is false an error will be raised when compound ids queried were not found in the database. 
		If true, just those that are found will be returned with no error or warning.
	}

}
\value{
   Returns a vector of compound names.The rownames will be the compound ids. Compound ids
   not found, or for which a name is not defined, will be represented as NA.
}
\author{
   Kevin Horan
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test2.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)

   getCompoundNames(conn,ids[1:3])
	unlink("test3.db")

}
