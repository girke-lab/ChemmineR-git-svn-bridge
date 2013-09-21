\name{smiles2sdf}
\alias{smiles2sdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SMILES (\code{character}) to \code{SDFset}
}
\description{
	Accepts a named vector or \code{SMIset} of 
	SMILES (Simplified Molecular Input Line Entry Specification) strings and
	returns its equivalent as an \code{SDFset} container. 
	
	This function runs in two modes. If ChemmineOB is available then it will use OpenBabel
	to convert all the given smiles into an SDFset with 2D coordinates. Otherwise
	the compound is submitted to the ChemMine Tools web service for conversion
	with the Open Babel Open Source Chemistry Toolbox. In this case only the first element
	will be used since this is a very slow operation.
}
\usage{
	smiles2sdf(smiles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smiles}{
	  A named vector of SMILES strings. The names will be used to name the SDF objects.
}
}
\value{
	\item{SDFset}{ for details see ?"SDFset-class"}
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
## convert to sdf
data(smisample)
(sdf <- smiles2sdf(smisample[1:4]))
}
}
\keyword{ utilities }
