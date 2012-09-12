\name{write.SDFsplit}
\alias{write.SDFsplit}
\title{
SDF split function
}
\description{
Splits SD Files into any number of smaller SD Files
}
\usage{
write.SDFsplit(x, filetag, nmol)
}
\arguments{
  \item{x}{
object of class \code{SDFset}, \code{SDFstr}
}
  \item{filetag}{
string to prepend to file names
}
  \item{nmol}{
integer specifying number of molecules in split SD files
}
}
\details{
To split an SD File into smaller ones, one can read the source file into R with \code{read.SDFstr}
and write out smaller ones with \code{write.SDFsplit}. Note: when importing big SD Files, \code{read.SDFstr} 
will be much faster than \code{read.SDFset}, and there is no need to go through an SDFset object instance
in this case.
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}
\seealso{
Import/export functions: \code{read.SDFset}, \code{read.SDFstr}, \code{read.SDFstr}, \code{read.SDFset}
}
\examples{
## Load sample data
library(ChemmineR)
data(sdfsample)

\dontrun{## Create sample SD File with 100 molecules
write.SDF(sdfsample, "test.sdf")

## Read in sample SD File 
sdfstr <- read.SDFstr("test.sdf")

## Run export on SDFstr object
write.SDFsplit(x=sdfstr, filetag="myfile", nmol=10)

## Run export on SDFset object
write.SDFsplit(x=sdfsample, filetag="myfile", nmol=10)

}
}
\keyword{ utilities }
