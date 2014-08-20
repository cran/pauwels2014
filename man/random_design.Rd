\name{random_design}
\alias{random_design}

\title{
	Simulates a randim design process.
}
\description{
	Simulates a randim design process.
}
\usage{
random_design(knobj, sample_function, exps, seed, credits = 5000, 
 file_to_save = NULL, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
  \item{sample_function}{
	A sample function that takes a knowledge list as argument and outputs a sample from the associated posterior in a design matrix.
}
  \item{exps}{
A dataframe representing possible experiments to be performed with their cost. See \code{\link{exps}}.
}
  \item{seed}{
	A random number generator seed.
}
  \item{credits}{
	Total credit to be spent.
}
  \item{file_to_save}{
	A file where the updated knowledge list should be saved at each step.
}
  \item{verbose}{
	Should the process print information about on going computation.
}
}
\details{
	This implements a strategy consisting in choosing experiments randomly. If a file name is provided, the resulting object will be saved at the corresponding location at each step of the process.
}
\value{
	An updated knowledge list.
}

\author{
Edouard Pauwels
}
\seealso{
\code{\link{knobjs}}, \code{\link{sample_function_multi_mod_weight}}, \code{\link{sample_function_single_mod}}, \code{\link{experiment_list1}}, \code{\link{observables}}, \code{\link{sample_function}}, \code{\link{exps}}
}
\examples{
data(exps)
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


## Decrease parameter values for the example
knobj$global_parameters$max_it <- 2
knobj$global_parameters$n_simu_weights <- 3
knobj$global_parameters$sample_burn_in <- 5
knobj$global_parameters$sample_to_keep1 <- 10
knobj$global_parameters$sample_to_keep2 <- 10
knobj$global_parameters$n_multi_mod_weight <- 2
knobj$global_parameters$final_sample <- 5
knobj$global_parameters$final_sample_design <- 5

## Run the random design (this takes quite some time)
#knobj <- random_design(knobj, sample_function_single_mod, exps, seed = 1, credits = 400)
}

\keyword{ Design strategies }
