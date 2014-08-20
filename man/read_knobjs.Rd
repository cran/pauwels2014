\name{read_knobjs}
\alias{read_knobjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Summarizes pre-computed results.
}
\description{
	Allows to reproduce the figures of the paper (see also vignette).
}
\usage{
read_knobjs(objs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{objs}{
	A table of strings representing names of objects given in \code{\link{knobjs}}.
	}
}
\details{
	The function computes the risk of the estimated posterior mean as a function of credit spent for each object.
}
\value{
A list of tables, one for each string. For each of them, the first line is the risk of the posterior sample mean and the second line is the credit spent.
}
\author{
Edouard Pauwels
}
\seealso{
  \code{\link{knobjs}}
}
\examples{
data(knobjs)
data(exps)
sapply(	
	1:length(knobjs),
	function(k){
		assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
	}
)

read_knobjs(paste("knobjActMult", 1:10, sep=""))
}
\keyword{ package }

