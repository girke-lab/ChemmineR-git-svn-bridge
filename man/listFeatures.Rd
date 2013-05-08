\name{listFeatures}
\alias{listFeatures}
\title{
	List Features
}
\description{
	List the available features in the given database. These features can
	be used in the \code{\link{findCompounds}} function.
}
\usage{
	listFeatures(conn)
}
\arguments{
  \item{conn}{
	  Database connection
	}
}
\value{
	A vector of character feature names.
}
\author{
	Kevin Horan
}


\seealso{
	\code{\link{findCompounds}} 
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test7.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample,fct=function(sdfset) cbind(mw=MW(sdfset)))
	listFeatures(conn) # produces c("mw")
	unlink("test7.db")
}
