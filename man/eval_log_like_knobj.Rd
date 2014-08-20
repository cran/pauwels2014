\name{eval_log_like_knobj}
\alias{eval_log_like_knobj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Posterior function.
}
\description{
	Computes the posterior value associated to a given parameter value for a given knowledge list.
}
\usage{
eval_log_like_knobj(theta, knobj, fail_incoming = F, fit = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
A parameter named numeric vector.
}
  \item{knobj}{
A knowledge list. See \code{\link{knobjs}}
}
  \item{fail_incoming}{
A boolean indicating wether an error message is given by the ode solver.
}
  \item{fit}{
A parameter to be passed to the likelihood function. It indicates wether further prior information about the smoothness of the dynamical time course should be considered. This is used to guide local search posterior maximization methods.
}
}
\details{
	The function computes the log prior first and then the likelihood associated to all the time course data found in the \code{knobj$datas} slot. The likelihood terms are summed. The prior term and the likelihood terms are weighted, weights being the inverse the number of observations they represent. This is necessary in order to give comparable contributions to low resolution and high resolution experiments.
}
\value{
	A numerical value if \code{fail_incoming == FALSE}. A list containing a \code{res} numerical slot and a \code{fail} boolean slot representing wether the ode solver failed or not.
}
\author{
Edouard Pauwels
}
\seealso{
\code{\link{log_prior}}, \code{\link{knobjs}}, \code{\link{eval_kn_log_like}}
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

eval_log_like_knobj(knobj$global_parameters$true_params_T, knobj)

}
\keyword{Simulation functions}
