\name{getCompoundFeatures}
\alias{getCompoundFeatures}
\title{
	Get Compound Features
}
\description{
	Get feature values for specific compounds.
}
\usage{
getCompoundFeatures(conn, compoundIds, featureNames, filename = NA, keepOrder = FALSE, allowMissing = FALSE, batchSize = 1e+05)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
	}
  \item{compoundIds}{
	  A vector of compound_id numbers from this database. These are not compound names. Features will be fetched for
	  each compound given here.
	}
  \item{featureNames}{
	  A vector of features to fetch the value for, for each given compound.
	}
  \item{filename}{
	  If given, dump the results into a comma seperated values (CSV) file instead of returning a data frame. This can
	  avoid some potential memory limits when fetching large sets of data.
	}
  \item{keepOrder}{
	  Ensure that the output order of values matches the order in which the compound ids where given. This will make things a
	  little slower, so should only be used where required.
	}
  \item{allowMissing}{
	  If false, raise an exception if a compound cannot be found, otherwise just silently ignore it and return data for
	  whatever compound were found.
	}
  \item{batchSize}{
	  The number of compounds to fetch in a single query. If you find your running out of memory you can try reducing this
	  values, as well as try writing the result to a file using the \code{filename} parameter.
	}
}
\value{
	If \code{filename} is not given, returns a data frame with the compound_id and any given feature names. Each row represents
	one compound. If \code{filename} is given a filename then no value is returned, but the given file is created.
}
\author{
	Kevin Horan
}
\examples{
	#create and initialize a new SQLite database
   conn = initDb("test1.db")

	data(sdfsample)

	#load data 
	ids=loadSdf(conn,sdfsample,
			  function(sdfset) 
					data.frame(MW = MW(sdfset),  rings(sdfset,type="count",upper=6, arom=TRUE))
			 )

	f = getCompoundFeatures(conn,ids,c("mw","rings"))

	unlink("test1.db")


}
