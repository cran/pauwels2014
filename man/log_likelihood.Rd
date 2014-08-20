\name{log_likelihood}
\alias{log_likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User defined likelihood function. 
}
\description{
Noise is assumed to be independent for each entry. The default likelihood assumes the same heteroscedastic noise as the model used in \code{\link{add_noise}}.
}
\usage{
log_likelihood(simu, simu_subset, data, fit = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simu}{
	Simulated time courses.
}
  \item{simu_subset}{
	Subset of the simulated time course which relates to observed data.
}
  \item{data}{
	Observed data.
}
  \item{fit}{
	Should smoothness prior information about simulation added to the prior or not?
}
}
\details{
The noise model is gaussian with variance of the form \code{(0.01 + 0.04 * m^2)} where \code{m} is the mean.
}
\value{
	A numerical value.
}
\author{
Edouard Pauwels
}
\examples{
data(experiment_list1)
data(observables)

## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

simu <- simulate_experiment(knobj$global_parameters$true_params_T, knobj, experiment_list1$nothing)

simu_subset <- simu[
 knobj$global_parameters$tspan \%in\% observables[["mrnaLow"]]$reso, 
 observables[["mrnaLow"]]$obs
]

data <- add_noise(simu_subset)

log_likelihood(simu, simu_subset, data) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model specific functions}
