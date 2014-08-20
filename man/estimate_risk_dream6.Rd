\name{estimate_risk_dream6}
\alias{estimate_risk_dream6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Expected risk estimation (comparison with litterature).
}
\description{
	Implements the strategy of Dream6 challenge winning team.
}
\usage{
estimate_risk_dream6(thetas, knobj, experiment_fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thetas}{
	A sample from the posterior associated to the knoweldge list \code{\link{knobjs}}.
}
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
  \item{experiment_fun}{
	A function that represents the molecular perturbation to be performed. See \code{\link{experiment_list1}}.
}
}
\details{
	The global variable \code{\link{observables}} should be defined.
}
\value{
A dataframe with the following columns
  \item{Measurement}{Factor representing possible measurements. See \code{\link{observables}}.}
  \item{Risk}{The risk associated to this measurement.}
  \item{Cost}{The cost associated to this measurement.}
}
\author{
	Edouard Pauwels
}
\examples{
data(experiment_list1)
data(observables)
data(knobjs)
sapply(	
	1:length(knobjs),
	function(k){
		assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
	}
)

knobjActMult1$global_parameters$n_simu_weights <- 3

estimate_risk_dream6(knobjActMult1$datas[[1]]$thetas[1:10,], 
	knobjActMult1, experiment_list1$nothing)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Risk estimation function }
% __ONLY ONE__ keyword per line
