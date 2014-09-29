\name{connections}
\alias{DUD}
\alias{DrugBank}
\title{
	Database Connections
}
\description{
	Get a connection to one of the pre-build compound databases.
}
\usage{
DUD()

DrugBank()
}
\value{
	A connection object to the ether the DUD or DrugBank database.
	This object must be passed to other functions which make use of the connection.
}
\author{
	Kevin Horan
}

\examples{
	\dontrun{
		dbConn = DrugBank()
	}
}
