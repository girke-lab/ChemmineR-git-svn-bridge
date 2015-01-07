\name{ExtSDF-class}
\Rdversion{1.1}
\docType{class}
\alias{ExtSDF-class}
\alias{getAtomAttr,ExtSDF-method}
\alias{getBondAttr,ExtSDF-method}
\alias{show,ExtSDF-method}

\title{Class "ExtSDF"}
\description{
	This is a subclass of \code{SDF} and thus inherits all the
	slots and methods from that class. It adds a list of extended attributes
	for atoms and bonds. These attributes can curretnly only be populated
	from a V3000 formatted SDF file.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExtSDF", ...)}. The
function \code{\link{read.SDFset}} will also return objects of this class if
the argument \code{extendedAttributes} is set to "TRUE".
}
\section{Slots}{
  \describe{
    \item{\code{extendedAtomAttributes}:}{Object of class \code{"list"} }
    \item{\code{extendedBondAttributes}:}{Object of class \code{"list"} }
  }
}
\section{Methods}{
  \describe{
    \item{getAtomAttr}{\code{signature(x = "ExtSDF")}: returns atom attributes as a named list of values}
    \item{getBondAttr}{\code{signature(x = "ExtSDF")}: returns bond attributes as a named list of values}
    \item{show}{\code{signature(object = "ExtSDF")}: prints summary of \code{SDF} as well as any defined extended
							 attributes for the atoms or bonds}
	 }
}
\references{
SDF V3000 format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
	Kevin Horan
}

\seealso{
Related classes: SDF, SDFset, SDFstr, AP, APset
}
\examples{
	showClass("ExtSDF")
}
\keyword{classes}
