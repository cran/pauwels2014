\name{simulate_experiment_no_transform}
\alias{simulate_experiment_no_transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Link to the ode solver.
}
\description{
Simulate the kinetics of the system using the ode solver of package \code{deSolve}.
}
\usage{
simulate_experiment_no_transform(theta, initial_conditions, knobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
	A parameter vector in physical space to be passed to the solver
}
  \item{initial_conditions}{
	A vector of initial conditions to be passed to the solver
}
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}
}
}
\details{
The solver parameters are given by the \code{global_parameters} slot of \code{knobj}. The solver will call an executable which is defined in the \code{src} directory of the package under the name \code{model0_simplified_mrna_rates}. The number and order of parameters in \code{theta} and initial conditions in \code{initial_conditions} should correspond to the definition of the dynamical system as described in this file. See the vignette and manual of package \code{deSolve}.
}
\value{
	A \code{deSolve} table with attributes.
}
\author{
Edouard Pauwels
}
\examples{
## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

temp <- simulate_experiment_no_transform(knobj$global_parameters$true_params, 
 knobj$global_parameters$initial_conditions, knobj)

head(temp)
}
\keyword{Simulation functions}
