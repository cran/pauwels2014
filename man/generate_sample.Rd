\name{generate_sample}
\alias{generate_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An implementation of the Metropolis Hasting algorithm
}
\description{
This is an implementation of MH algorithm to sample from the posterior distribution. The proposal is a mixture between a gaussian proposal and a single coordinate proposal. The step size is diminished when the rejection rate is too high.
}
\usage{
generate_sample(theta, knobj, N = 500, step = 1, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  An initialization parameter named numeric vector
}
  \item{knobj}{
A knowledge list. See \code{\link{knobjs}}.
}
  \item{N}{
The total sample size.
}
  \item{step}{
The proposal distribution expected step length.
}
  \item{verbose}{
Should the sampling process print information about itself?
}
}
\details{
The posterior is evaluated using \code{\link{eval_log_like_knobj}} function.
}
\value{
A posterior sample matrix, each row representing a parameter named numeric vector.
}
\author{
Edouard Pauwels
}

\seealso{
\code{\link{eval_log_like_knobj}}, \code{\link{knobjs}}
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

generate_sample(knobj$global_parameters$params * 50, knobj, N = 10)

}
\keyword{Sampling functions}

