\name{sdf2smiles}
\alias{sdf2smiles}
\alias{sdf2smilesOB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{SDFset} to \code{character}
Convert \code{SDFset} to SMILES (\code{character})
}
\description{
Accepts compounds in an \code{SDFset} container and returns the corresponding
SMILES (Simplified Molecular Input Line Entry Specification) strings as \code{SMIset} object. 
If ChemineOB is available then OpenBabel for the format conversion. 
Otherwise the compound is submitted to the ChemMine Tools web service for conversion
with the Open Babel Open Source Chemistry Toolbox. If the input object
contains multiple items, only the first is converted.
}
\usage{
	sdf2smiles(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
  	  A \code{SDFset} object which containing the given compounds
}
}
\value{
	\item{character}{ for details see ?"character"}
}
\references{
Chemmine web service: http://chemmine.ucr.edu

Open Babel: http://openbabel.org

SMILES Format: http://en.wikipedia.org/wiki/Chemical_file_format#SMILES
}
\author{
Tyler Backman, Kevin Horan
}

\examples{
\dontrun{
## get a sample compound
data(sdfsample); sdfset <- sdfsample[1]
## convert to smiles
(smiles <- sdf2smiles(sdfset))
as.character(smiles)
}
}
\keyword{ utilities }
