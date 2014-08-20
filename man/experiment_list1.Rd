\name{experiment_list1}
\alias{experiment_list1}
\docType{data}
\title{
	Molecular perturbations.
}
\description{
	A list of functions that represent possible molecular perturbations to be performed on the network. 
}
\usage{data(experiment_list1)}
\details{
A list of 10 functions. Each function takes as input a named numeric parameter vector \code{theta} and an initial condition named numeric vector \code{initial_conditions}. Each function outputs a list with the following slots:
\itemize{
	\item \code{theta}: A modified parameter value
	\item \code{initial_conditions}: A modified initial condition value
	\item \code{cost}: The cost of the associated molecular perturbation
}
The 10 functions perform the following operations on their input (parameters numbering might not be always coherent with names):
\itemize{
	\item \code{delete_gene6}
	\itemize{
		\item \code{theta}: set \code{pro6_strength} and \code{rbs6_strength}
		\item \code{initial_condtions}: set \code{p6} initial concentration to 0
		\item \code{cost}: 800
	}
	\item {delete_gene7} 
	\itemize{
		\item \code{theta}: set \code{pro7_strength} and \code{rbs8_strength}
		\item \code{initial_condtions}: set \code{p7} initial concentration to 0
		\item \code{cost}: 800
	}
	\item \code{delete_gene8} 
	\itemize{
		\item \code{theta}: set \code{pro7_strength} and \code{rbs8_strength}
		\item \code{initial_condtions}: set \code{p7} initial concentration to 0
		\item \code{cost}: 800
	}
	\item \code{knockdown_gene6} 
	\itemize{
		\item \code{theta}: multiply \code{mrna6_degradation_rate} by 10
		\item \code{initial_condtions}: Nn operation
		\item \code{cost}: 350
	}
	\item \code{knockdown_gene7} 
	\itemize{
		\item \code{theta}: multiply \code{mrna7_degradation_rate} by 10
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 350
	}
	\item \code{knockdown_gene8} 
	\itemize{
		\item \code{theta}: multiply \code{mrna8_degradation_rate} by 10
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 350
	}
	\item \code{decrease_rbs_gene6} 
	\itemize{
		\item \code{theta}: divide \code{rbs6_strength} by 10
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 350
	}
	\item \code{decrease_rbs_gene7} 
	\itemize{
		\item \code{theta}: divide \code{rbs6_strength} by 10
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 350
	}
	\item \code{decrease_rbs_gene8} 
	\itemize{
		\item \code{theta}: divide \code{rbs7_strength} by 10
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 350
	}
	\item \code{nothing}
	\itemize{
		\item \code{theta}: no operation
		\item \code{initial_condtions}: no operation
		\item \code{cost}: 0
	}
}
}
\examples{
data(experiment_list1)
data(knobjs)
sapply(	
	1:length(knobjs),
	function(k){
		assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
	}
)

theta <- knobjActMult1$datas[[1]]$thetas[1,]
thetaT <- knobjActMult1$transform_params(theta)

temp <- experiment_list1$delete_gene7(thetaT, knobjActMult1$global_parameters$initial_conditions)

rbind(temp$theta, thetaT)
rbind(temp$initial_conditions, knobjActMult1$global_parameters$initial_conditions)

}
\keyword{datasets}
\keyword{Model specific functions}
