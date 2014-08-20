\name{observables}
\alias{observables}
\docType{data}
\title{
Observable quantities of the model
}
\description{
A list of quantities that can be observed. Each entry describes the quantity to observe and its time resolution as well as the cost of the corresponding observation.
}
\usage{data(observables)}
\format{
A list with the follwing entries: \code{"p6"}, \code{"p7"}, \code{"p8"}, \code{"mrnaHigh"}, \code{"mrnaLow"}. Each of those entries is a list with the following slots:
\itemize{
	\item \code{name}: (chr) name of the observable
	\item \code{obs}: (chr) the subset of kinetic variables to observe
	\item \code{reso}: (num) time points to observe
	\item \code{cost}: (nul) the cost of observing this quantity
}
}
\details{
The complete description of the list is as follows
\itemize{
	\item \code{"p6"}:
	\itemize{
		\item \code{name}: \code{"p6"}
		\item \code{obs}: \code{c("time", "p6")}
		\item \code{reso}: from \code{0} to \code{100} by steps of \code{0.5}
		\item \code{cost}: 400
	}
	\item \code{"p7"}:
	\itemize{
		\item \code{name}: \code{"p7"}
		\item \code{obs}: \code{c("time", "p7")}
		\item \code{reso}: from \code{0} to \code{100} by steps of \code{0.5}
		\item \code{cost}: 400
	}
	\item \code{"p8"}:
	\itemize{
		\item \code{name}: \code{"p8"}
		\item \code{obs}: \code{c("time", "p8")}
		\item \code{reso}: from \code{0} to \code{100} by steps of \code{0.5}
		\item \code{cost}: 400
	}
	\item \code{"mrnaHigh"}:
	\itemize{
		\item \code{name}: \code{"mrnaHigh"}
		\item \code{obs}: \code{c("time", "v6_mrna", "v7_mrna", "v8_mrna")}
		\item \code{reso}: from \code{0} to \code{100} by steps of \code{2}
		\item \code{cost}: 1000
	}
	\item \code{"mrnaLow"}:
	\itemize{
		\item \code{name}: \code{"mrnaLow"}
		\item \code{obs}: \code{c("time", "v6_mrna", "v7_mrna", "v8_mrna")}
		\item \code{reso}: from \code{0} to \code{100} by steps of \code{4}
		\item \code{cost}: 500
	}
}
}
\examples{
data(observables)
data(knobjs)
sapply(	
	1:length(knobjs),
	function(k){
		assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
	}
)
data(exps)

theta <- knobjActMult1$datas[[1]]$thetas[1,]
thetaT <- knobjActMult1$transform_params(theta)

temp <- simulate_experiment_no_transform(thetaT, 
	knobjActMult1$global_parameters$initial_conditions, knobjActMult1)

observable <- observables$mrnaLow
temp[temp[,1] \%in\% observable$reso, colnames(temp) \%in\% observable$obs]

}
\keyword{datasets}
