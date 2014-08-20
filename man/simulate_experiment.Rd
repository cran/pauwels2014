\name{simulate_experiment}
\alias{simulate_experiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates the dynamics of a molecular perturbation
}
\description{
This function simulates the kinetics of the system, in accordance with a chosen molecular perturbation.
}
\usage{
simulate_experiment(theta, knobj, experiment_fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
A named numeric parameter vector
}
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
  \item{experiment_fun}{
	An element of \code{\link{experiment_list1}} that represents a molecular perturbation
}
}
\details{
Transforms parameter values and initial conditions in order to pass it to \code{deSolve} ode solver.
}
\author{
Edouard Pauwels
}
\seealso{
\code{\link{transform_params}}, \code{\link{experiment_list1}}, \code{\link{simulate_experiment_no_transform}}
}
\examples{
data(experiment_list1)

## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

temp <- simulate_experiment(knobj$global_parameters$true_params, 
 knobj, experiment_list1$nothing)

head(temp)
}
\keyword{Simulation functions}
