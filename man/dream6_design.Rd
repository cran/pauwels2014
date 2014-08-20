\name{dream6_design}
\alias{dream6_design}
\title{
	Simulates the active design process using the comparison criterion (see article for details).
}
\description{
	Simulates the active design process.
}
\usage{
dream6_design(knobj, sample_function, 
	seed, credits = 5000, file_to_save = NULL, verbose = T)
}
\arguments{
  \item{knobj}{
	A knowledge list. See \code{\link{knobjs}}.
}
  \item{sample_function}{
	A sample function that takes a knowledge list as argument and outputs a sample from the associated posterior in a design matrix.
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
	This implements the active design procedure using the criterion used by the wining strategy for DREAM6 challenge. If a file name is provided, the resulting object will be saved at the corresponding location at each step of the process. The function requires the global variables \code{\link{experiment_list1}} and \code{\link{observables}} to be available.
}
\value{
	An updated knowledge list.
}
\author{
	Edouard Pauwels
}
\seealso{
	\code{\link{knobjs}}, \code{\link{sample_function_multi_mod_weight}}, \code{\link{sample_function_single_mod}}, \code{\link{experiment_list1}}, \code{\link{observables}}, \code{\link{sample_function}}.
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
knobj$global_parameters$n_simu_weights <- 2
knobj$global_parameters$sample_burn_in <- 5
knobj$global_parameters$sample_to_keep1 <- 2
knobj$global_parameters$n_multi_mod <- 2
knobj$global_parameters$final_sample <- 2
knobj$global_parameters$final_sample_design <- 2

## Run the active design (this takes quite some time)
#knobj <- dream6_design(knobj, sample_function_single_mod, seed = 1, credits = 400)

}

\keyword{ Design strategies }
