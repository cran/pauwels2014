\name{eval_kn_log_like}
\alias{eval_kn_log_like}
\title{
Evaluates the likelihood of a parameter value
}
\description{
This evaluates the likelihood of a parameter value by comparing corresponding kinetics to a given time course matrix. 
}
\usage{
eval_kn_log_like(theta, initial_conditions, data, knobj, 
 fail_incoming = F, simu = NULL, fit = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
A parameter named numeric vector.
}
  \item{initial_conditions}{
Initial conditions named numeric vector.
}
  \item{data}{
A time course matrix
}
  \item{knobj}{
A knowledge list. See \code{\link{knobjs}}
}
  \item{fail_incoming}{
A boolean indicating wether an error message is given by the ode solver.
}
  \item{simu}{
The simulated time course matrix corresponding to the parameter \code{theta}. If it is not provided, it will be computed by a call to \code{\link{simulate_experiment_no_transform}}.
}
  \item{fit}{
A parameter to be passed to the likelihood function. It indicates wether further prior information about the smoothness of the dynamical time course should be considered. This is used to guide local search posterior maximization methods.
}
}
\details{
The comparison is made based on the user defined \code{\link{log_likelihood}} function.
}
\value{
A numerical value if \code{fail_incoming == FALSE}. A list containing a \code{res} numerical slot and a \code{fail} boolean slot representing weither the ode solver failed or not.
}
\author{
Edouard Pauwels
}
\seealso{
\code{\link{log_likelihood}}, \code{\link{knobjs}}, \code{\link{simulate_experiment_no_transform}}
}
\examples{
data(experiment_list1)
data(observables)

## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

## Initialize with some data
knobj$datas[[1]] <- list(
 manip = experiment_list1$nothing,
 data = add_noise(
  simulate_experiment(knobj$global_parameters$true_params_T, knobj, experiment_list1$nothing)[
   knobj$global_parameters$tspan \%in\% observables[["mrnaLow"]]$reso, 
   observables[["mrnaLow"]]$obs
  ]
 )
)

eval_kn_log_like(
 knobj$global_parameters$true_params, 
 knobj$global_parameters$initial_conditions, 
 knobj$datas[[1]]$data, knobj )

}

\keyword{Simulation functions}
