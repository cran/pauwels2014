\name{risk_theta_vect}
\alias{risk_theta_vect}

\title{
Expected risk based on a posterior sample
}
\description{
The function computes the expected risk based on the empirical distribution defined by the sample thetas_trans.
}
\usage{
risk_theta_vect(thetas_trans, n_params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thetas_trans}{
A matrix which row represent parameter values in physical space (after applying \code{\link{transform_params}}).
}
  \item{n_params}{
The number of parameters to be estimated
}
}
\value{
A numerical value
}
\author{
Edouard Pauwels
}
\seealso{
\code{\link{risk_theta_fun}}
}
\examples{
data(knobjs)
sapply(	
	1:length(knobjs),
	function(k){
		assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
	}
)

thetas <- knobjActMult1$datas[[1]]$thetas[1:10,]

thetas_trans <- t(apply(thetas, 1, transform_params))
risk_theta_vect(thetas_trans, 9)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model specific functions}
\keyword{Risk estimation functions}
