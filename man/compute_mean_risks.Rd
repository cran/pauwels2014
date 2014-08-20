\name{compute_mean_risks}
\alias{compute_mean_risks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute an average risk as a function of credit spent
}
\description{
	Used to summarize information contained in lists given by \code{\link{read_knobjs}}.
}
\usage{
compute_mean_risks(mean_risks, legend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean_risks}{
		A list of tables as given by \code{\link{read_knobjs}}.
}
  \item{legend}{
	A label to be associated to the output data frame.
}
}
\details{
	Summarizes results in a data frame for further ploting.
}
\value{
	A data frame with four columns: \code{cost}, the budget spent, numerical, \code{risk}, the associated risk, numerical, \code{type}, given by the \code{legend} argument and \code{chain} which is an integer associated to each table in the list argument.
}
\author{
Edouard Pauwels
}
\seealso{
 \code{\link{read_knobjs}} 
}
\examples{
	data(knobjs)
	sapply(	
		1:length(knobjs),
		function(k){
			assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
		}
	)

	data(exps)
	temp <- read_knobjs(paste("knobjActMult", 1:10, sep=""))
	mean_risks <- compute_mean_risks(temp, "A title")
	mean_risks

}

\keyword{ Optimization functions }

