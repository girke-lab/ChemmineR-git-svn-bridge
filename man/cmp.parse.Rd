\name{cmp.parse}
\alias{cmp.parse}
\title{Parse an SDF file and compute descriptors for all compounds}
\description{
'cmp.parse' will take a SDF file, parse all the compounds
encoded, compute their atom-pair descriptors, and return the descriptors as a
list. The list contains two names, 'descdb' and 'cids'. 'descdb' is a vector of
descriptors, and 'cids' is a list of names of compounds found in the SDF file.
The returned list is usually used to a database, against which similarity
search can be performed using the 'search' function. These two functions will
parse all compounds in the SDF file. To parse a single compound, use
'cmp.parse1' instead.
}
\usage{
cmp.parse(filename, quiet=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{The file name of the SDF file}
  \item{quiet}{Whether to silent the output of progress information}
}
\details{
	The 'filename' can be a local file or an URL. It is interactive, and
	will display the parsing progress. Since the parsing will also compute
	of atom-pair descriptors, it is time consuming. You will be reminded to
	save the parsing result for future use at the end of parsing.
}
\value{
	Return a list that can be used as the database against which similarity
	search can be performed. The 'search' and 'cmp.cluster' functions both
	expect a database returned by 'cmp.parse'.
  \item{descdb}{A vector containing the descriptors for all the compounds.}
  \item{cids}{Compound ID information found in the SDF file. It is the first
	  line of SDF of a compound.}
}
\references{Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", in \emph{J Chem Inf Comput Sci}.}
\author{Y. Eddie Cao, Li-Chang Cheng}

\seealso{\code{\link{cmp.parse1}}, \code{\link{cmp.search}},
	\code{\link{cmp.cluster}},
\code{\link{cmp.similarity}}}
\examples{
# load sample database from web
db <- cmp.parse("http://bioweb.ucr.edu/ChemMineV2/static/example_db.sdf")
# (optinally) save the db for future use
save(db, file="db.rda", compress=TRUE)
# ...
# later, in a separate session, you can load it back:
load("db.rda")
}
\keyword{utilities}
