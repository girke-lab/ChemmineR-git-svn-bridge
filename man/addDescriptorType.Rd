\name{addDescriptorType}
\alias{addDescriptorType}
\title{
	Add Descriptor Type
}
\description{
	Add a new descriptor type to the database. Normally descriptor types are added as needed, but
	if you are doing a parrallel data load you must pre-load the descriptor type to prevent
	duplicate defintion errors.
}
\usage{
addDescriptorType(conn, descriptorType)
}
\arguments{
  \item{conn}{
		Any database connection object.
	}
  \item{descriptorType}{
	  The name of the descriptor.
	}
}
\value{
	No return value.
}
\author{
	Kevin Horan
}
\examples{
	\dontrun{
		conn =  initDb(...)
		addDescriptor(conn,"fp")
	}
}
