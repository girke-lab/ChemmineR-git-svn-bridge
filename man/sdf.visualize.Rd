\name{sdf.visualize}
\alias{sdf.visualize}
\title{Visualize an SDFset online using ChemMine Tools}
\description{
'sdf.visualize' will take an \code{SDFset} object and send the compounds to the ChemMine Tools
website, for visualization and futher analysis. The results are launched in the users web browser.
}
\usage{
sdf.visualize(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{sdf}{
	  A \code{SDFset} object which containing the given compounds
	}
}
\references{
ChemMine Tools web service: http://chemmine.ucr.edu
}
\author{
Tyler Backman
}
\value{
	Returns the URL of the webpage containing all the SDFs and 2D images
	corresponding to the selected compounds.
}

\seealso{\code{\link{cmp.parse}}, \code{\link{sdf.subset}}, \code{\link{plotStruc}}}
\examples{
## Load sample SD file
data(sdfsample)
sdfset <- sdfsample

\dontrun{
## Plot structures using web service ChemMine Tools
sdf.visualize(sdfset[1:4])
}
}
\keyword{utilities}
