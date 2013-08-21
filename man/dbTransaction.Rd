\name{dbTransaction}
\alias{dbTransaction}
\title{
	DB Transaction
}
\description{
	Run any db statements inside a transaction. If any error is raised the transaction will
	be rolled back, otherwise it will be committed at the end. 
}
\usage{
dbTransaction(conn, expr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
	}

  \item{expr}{
		Any block of code.
	}
}
\value{
	The value of the given block of code will be returned upon successfully commiting the	
	transaction. Otherwise an error will be raised.
}
\author{
	Kevin Horan
}
\examples{
	
   conn = initDb("test15.db")
	dbTransaction(conn,{
		# any db code here
	})
}
