\name{sample_function_multi_mod_weight}
\alias{sample_function_multi_mod_weight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample function visiting multiple modes of the posterior
}
\description{
Generate a posterior sample using multiple local search maximization and sampling based on different initializations.
}
\usage{
sample_function_multi_mod_weight(knobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
}
\details{
	The parameters governing the local search and sampling behaviour are defined in the \code{global_parameters} slot of the \code{\link{knobjs}} argument. The function consists in using the \code{\link{BFGS_special}} function to find an initialization for the Metropolis Hasting algorithm implented by \code{\link{generate_sample}}. This is done multiple times. This procedure is also applied to previous sample points. All those sample are aggregated and the resulting sample is chosen randomly based on the associated posterior values.
}
\value{
A matrix which rows represent a named numeric vector of parameters
}
\author{
Edouard Pauwels
}
\seealso{
\code{\link{sample_function}}, \code{\link{BFGS_special}}, \code{\link{generate_sample}}
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

## Decrease parameter values for the example
knobj$global_parameters$max_it <- 2
knobj$global_parameters$n_multi_mod <- 2
knobj$global_parameters$sample_burn_in <- 5
knobj$global_parameters$sample_to_keep1 <- 100
knobj$global_parameters$final_sample <- 100
knobj$global_parameters$final_sample_est <- 100


#thetas <- sample_function_multi_mod_weight(knobj)
#thetas
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sampling functions}

