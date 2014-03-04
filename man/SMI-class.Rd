\name{SMI-class}
\Rdversion{1.1}
\docType{class}
\alias{SMI-class}
\alias{as.character,SMI-method}
\alias{coerce,character,SMI-method}
\alias{show,SMI-method}

\title{Class \code{"SMI"}}
\description{
Container for storing the SMILES string of a single molecule.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SMI", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{smiles}:}{Object of class \code{"character"} of length one}
  }
}
\section{Methods}{
  \describe{
    \item{as.character}{\code{signature(x = "SMI")}: returns content as character vector}
    \item{coerce}{\code{signature(from = "character", to = "SMI")}: \code{as(smi, "SMI")}}
    \item{coerce}{\code{signature(from = "SMIset", to = "SMI")}: \code{as(smiset, "SMI")}}
    \item{show}{\code{signature(object = "SMI")}: prints summary of \code{SMI} }
	 }
}
\references{
SMILES (Simplified molecular-input line-entry system) format definition:
http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
}
\author{
Thomas Girke
}
\seealso{
Related classes: SMIset, SDF, SDFset
}
\examples{
showClass("SMI")

## Instances of SMI class
data(smisample); smiset <- smisample
(smi <- smiset[[1]]) # returns first molecule in smiset as SMI object
}
\keyword{classes}
