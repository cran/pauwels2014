\name{exps}
\alias{exps}
\docType{data}
\title{
	List of possible experiments
}
\description{
	An experiment associates a molecular perturbation and a quantity to observe.
}
\usage{data(exps)}
\format{
  A data frame that combines all possible experiments of \code{\link{experiment_list1}} to all possible observations of \code{\link{observables}}.
  }
\details{
The data frame contains 49 entries of the following variables
\itemize{
    \item \code{Measurement}: a factor with levels \code{"mrnaHigh"}, \code{"mrnaLow"}, \code{"p6"}, \code{"p7"}, \code{"p8"}. See \code{\link{observables}}
    \item \code{Cost}: a numeric vector. The cost of the corresponding experiment and observation to be performed.
    \item \code{exp}: a factor with levels \code{"decrease_rbs_gene6"}, \code{"decrease_rbs_gene7"}, \code{"decrease_rbs_gene8"}, \code{"delete_gene6"}, \code{"delete_gene7"}, \code{"delete_gene8"}, \code{"knockdown_gene6"}, \code{"knockdown_gene7"}, \code{"knockdown_gene8"}, \code{"nothing"}. See \code{\link{experiment_list1}} for a description of the corresponding perturbations
  }
This should be 50 experiments by default (10 experiments and 5 observables). However, we always initialise the design process with experiment \code{nothing} and observable \code{mrnaLow}, this leaves 49 possible experiments left. This object is used by the \code{\link{random_design}} function. 
}
\examples{
data(exps)
exps[1,]
}
\keyword{datasets}
\keyword{Model specific objects}
