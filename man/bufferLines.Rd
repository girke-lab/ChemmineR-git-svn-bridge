\name{bufferLines}
\alias{bufferLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Buffer File Input
}
\description{
	Buffer the input of files to increase efficiency
}
\usage{
bufferLines(fh, batchSize, lineProcessor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fh}{
	  file handle
	}
  \item{batchSize}{
		How many lines to read in each batch
	}
  \item{lineProcessor}{
		Each batch of lines will be passed to this function for processing
	}
}
\value{
	No return value
}
\author{
	Kevin Horan
}

\examples{

	\dontrun{
		fh = file("filename")
		bufferLines(fh,100,function(lines) {
			message("found ",length(lines)," lines")
		})
	}
}
