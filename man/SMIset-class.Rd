\name{SMIset-class}
\Rdversion{1.1}
\docType{class}
\alias{SMIset-class}
\alias{[,SMIset-method}
\alias{[[,SMIset-method}
\alias{[<-,SMIset-method}
\alias{as.character,SMIset-method}
\alias{c,SMIset-method}
\alias{cid,SMIset-method}
\alias{cid<-,SMIset-method}
\alias{coerce,character,SMIset-method}
\alias{coerce,list,SMIset-method}
\alias{coerce,SMIset,SMI-method}
\alias{length,SMIset-method}
\alias{show,SMIset-method}
\alias{view,SMIset-method}

\title{Class \code{"SMIset"}}
\description{
List-like container for storing SMILES strings of many compounds. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SMIset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{smilist}:}{Object of class \code{"list"} with compound identifiers stored in name slots}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SMIset")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SMIset")}: returns single component as \code{SMI} object }
    \item{[<-}{\code{signature(x = "SMIset")}: replacement method for one or many entries }
    \item{as.character}{\code{signature(x = "SMIset")}: returns content as named character vector }
    \item{c}{\code{signature(x = "SMIset")}: concatenates two \code{SMIset} containers }
    \item{cid}{\code{signature(x = "SMIset")}: returns compound identifiers }
    \item{cid<-}{\code{signature(x = "SMIset")}: replacement method for compound identifiers }
    \item{coerce}{\code{signature(from = "character", to = "SMIset")}: \code{as(character, "SMIset")} }
    \item{coerce}{\code{signature(from = "list", to = "SMIset")}: \code{as(list, "SMIset")} }
    \item{coerce}{\code{signature(from = "SMIset", to = "SMI")}: \code{as(smiset, "SMI")} }
    \item{length}{\code{signature(x = "SMIset")}: returns number of entries stored in object }
    \item{show}{\code{signature(object = "SMIset")}: prints summary of \code{SMIset} }
    \item{view}{\code{signature(x = "SMIset")}: prints extended summary of \code{SMIset} }
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
Related classes: SMI, SDF, SDFset

Import function: read.SMIset("some_SMILES_file") 

Export function: write.SMI(smiset, "some_file.smi")
}
\examples{
showClass("SMIset")

## Instances of SMIset class
data(smisample); smiset <- smisample
smiset; view(smiset[1:4])
smiset[[1]]

## Import and store SMILES file in SMIset container
# smiset <- read.SMIset("some_SMILES_file")

## Miscellaneous accessor methods
cid(smiset[1:4])
(smivec <- as.character(smiset[1:4]))

## Construct SMIset from named vector 
as(smivec, "SMIset")

## Assigning compound IDs and keeping them unique
unique_ids <- makeUnique(cid(smiset))
cid(smiset) <- unique_ids

## Export SMIset to SMILES file
# write.SMI(smiset[1:4], file="sub.smi", cid=TRUE)

}
\keyword{classes}
