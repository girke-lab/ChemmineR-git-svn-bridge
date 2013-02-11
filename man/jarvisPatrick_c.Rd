\name{jarvisPatrick_c}
\alias{jarvisPatrick_c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Jarvis Patrick Clustering in C code
}
\description{
	This not meant to be used directly, use jarvisPatrick instead. It is exposed so other
	libraries can make use of it.
}
\usage{
jarvisPatrick_c(neighbors, minNbrs, fast = TRUE, bothDirections = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{neighbors}{
	  A matrix of integers. Non integer matricies will be coerced. 
	  Each row represensts one element, indexed 1 to N. The values in row
	  i should be the index value of the neighbors of i. Thus, each value should 
	  itself be a valid row index. 
}
  \item{minNbrs}{
	  The minimum number of common neibhbors needed for two elements to be merged.
}
  \item{fast}{
	  If true, only the neibhors given in each row are checked to see if they share \code{minNbrs} 
	  neighbors in common. If false, all pairs of elements are compared. For a matrix of size NxM,
	  the first method yeilds a running time of O(NM), while the second yeilds a running time of O(N^2).
}
  \item{bothDirections}{
	  If true, two elements must contain each other in their neighbor list in order to be merged.
	  If false and fast is true, then only one element must contain the other as a neighbor. If
	  false and fast is false, than neither element must contain the other as a neighbor, though
	  in all cases there must still be at least \code{minNbrs} neibhros in common.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (neighbors, minNbrs, fast = TRUE, bothDirections = FALSE) 
{
    n = matrix(as.integer(neighbors), nrow(neighbors), ncol(neighbors))
    .Call("jarvis_patrick", n, as.integer(minNbrs), as.integer(fast), 
        as.integer(bothDirections))
}
}
