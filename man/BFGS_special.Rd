\name{BFGS_special}
\alias{BFGS_special}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	An implementation of BFGS method for posterior maximization.
}
\description{
	Gradients are computed using finite differences.
}
\usage{
BFGS_special(init, knobj, fun_like, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
	An initial value of the parameter to be optimized.
}
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
  \item{fun_like}{
	A function to compute posterior value. See \code{\link{eval_log_like_knobj}}
}
  \item{verbose}{
	Print progresses of the local search?
}
}
\details{
	The step size are chosen using Armijo's rule. Special checks are performed to avoid numerical instabilities in the differential equation solver. 
}
\value{
A list with the following entries:
  \item{theta}{The local optimum found by the method.}
  \item{fail}{A boolean representing wither the local search failed or not due to numerical problems.}
%% ...
}
\author{
Edouard Pauwels
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
knobj$experiments <- paste("nothing", "mrnaLow")

theta <- rep( 50, length(knobj$global_parameters$param_names) )
names(theta) <- knobj$global_parameters$param_names

## Only perform 5 iterations
knobj$global_parameters$max_it <- 5

temp <- BFGS_special(theta, knobj, eval_log_like_knobj)
temp$theta
}

\keyword{ Optimization functions }

