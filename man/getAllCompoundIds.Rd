\name{getAllCompoundIds}
\alias{getAllCompoundIds}
\title{
	Get ALl Compound Ids
}
\description{
	Return a vector of every compound id in the given database. 
}
\usage{
getAllCompoundIds(conn)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
	}
}
\value{
	A vector of compound_id numbers
}
\author{
	Kevin Horan
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test1.db")

	data(sdfsample)

	#load data
	ids=loadSdf(conn,sdfsample)
	ids2=getAllCompoundIds(conn)
	#ids == ids2

	unlink("test1.db")
}
