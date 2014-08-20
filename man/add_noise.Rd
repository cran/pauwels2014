\name{add_noise}
\alias{add_noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Noise generative process for the simulations
}
\description{
	Specifies a noise generative process for the simulations. This describes how the true dynamics of the system is perturbed by noise.	
}
\usage{
add_noise(data_theta_Ts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_theta_Ts}{
	A time series data matrix, the first column representing time and the remaining columns representing time course of various quantities of interest.
}
}
\details{
	The default generative process is to add independant gaussian heteoscedastic noise to all columns, except the first one representing time. The noise model is gaussian with variance of the form \code{(0.01 + 0.04 * m^2)} where \code{m} is the mean.
}
\value{
	A time series data matrix of the same size as the input.
}
\author{
	Edouard Pauwels
}
\examples{
data(experiment_list1)
data(observables)

## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

## Generate a time cours matrix
tempCourse <- simulate_experiment(
 knobj$global_parameters$true_params_T, 
 knobj, 
 experiment_list1$nothing
)[
 knobj$global_parameters$tspan \%in\% observables[["mrnaLow"]]$reso, 
 observables[["mrnaLow"]]$obs
]

## Add noise to the time course matrix
add_noise(tempCourse)
}
\keyword{Model specific functions}
