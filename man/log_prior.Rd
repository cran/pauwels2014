\name{log_prior}
\alias{log_prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User defined log prior
}
\description{
	User defined log prior based on the parameter value only.
}
\usage{
log_prior(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
	A parameter named numeric vector.
}
}
\details{
	The default is a  gaussian prior with large covariance matrix. This is actually a prior on the parameters in log space. The prior can can thus be thought of as a log normal prior.
}
\value{
A numerical value
}
\author{
Edouard Pauwels
}
\examples{
theta <- 50 + sample(c(0,1), size = 34, replace = TRUE)
log_prior(theta)
}
\keyword{Model specific functions}
