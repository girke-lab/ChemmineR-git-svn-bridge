\name{FP-class}
\Rdversion{1.1}
\docType{class}
\alias{FP-class}
\alias{as.character,FP-method}
\alias{as.numeric,FP-method}
\alias{as.vector,FP-method}
%\alias{coerce,FPset,FP-method}
\alias{coerce,numeric,FP-method}
\alias{show,FP-method}
\alias{c,FP-method}
\alias{fold,FP-method}
\alias{foldCount,FP-method}
\alias{fptype,FP-method}
\alias{numBits,FP-method}

\title{Class \code{"FP"}}
\description{
Container for storing the fingerprint of a single compound. The \code{FPset} class is used for storing the fingerprints of many compounds.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FP", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fp}:}{Object of class \code{"numeric"} }
    \item{\code{foldCount}:}{Object of class \code{"numeric"} }
    \item{\code{type}:}{Object of class \code{"character"} }
  }
}
\section{Methods}{
  \describe{
    \item{as.character}{\code{signature(x = "FP")}: returns fingerprint as character string  }
    \item{as.numeric}{\code{signature(x = "FP")}: returns fingerprint as numeric vector }
    \item{as.vector}{\code{signature(x = "FP")}: returns fingerprint as numeric vector }
    \item{coerce}{\code{signature(from = "FPset", to = "FP")}: coerce \code{FPset} component to list with many \code{FP} objects }
    \item{coerce}{\code{signature(from = "numeric", to = "FP")}: construct \code{FP} object from numeric vector }
    \item{show}{\code{signature(object = "FP")}: prints summary of \code{FP} }
    \item{c}{\code{signature(x = "FP")}: concatenates any number of \code{FP} objects}
    \item{fold}{\code{signature(x = "FP")}: fold fingerprint in half }
    \item{foldCount}{\code{signature(x = "FP")}: number of times this object has been folded}
    \item{fptype}{\code{signature(x = "FP")}: the type of this fingerprint }
    \item{numBits}{\code{signature(x = "FP")}: the number of bits in this fingerprint}
	 }
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", in J Chem Inf Comput Sci.
}
\author{
Thomas Girke
}
\seealso{
Related classes: SDF, SDFset, SDFstr, AP, APset, FPset. 
}
\examples{
showClass("FP")

## Instance of FP class
data(apset)
fpset <- desc2fp(apset)
(fp <- fpset[[1]])

## Class usage
fpc <- as.character(fp)
fpn <- as.numeric(fp)
as(fpn, "FP")
as(fpset[1:4], "FP") 

}
\keyword{classes}
