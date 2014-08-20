\name{log_normalize}
\alias{log_normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Normalize in log space
}
\description{
	Binding to a \code{C} implementation of addition and normalization in log space.
}
\usage{
log_normalize(vals, maxi = max(vals), ID = which.max(vals)[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vals}{
	Log-values to be normalize 
}
  \item{maxi}{
	The maximum value of \code{vals}
}
  \item{ID}{
	One entry of \code{vals} which has the highest value.
}
}
\value{
	A list with three slots. First one is a numerical vector with values normalized so that they sum to 1. The second is the same normalized vector in log space. The last one is the log of the normalization constant. See example.	
}

\author{
Edouard Pauwels
}
\examples{
log_x <- rnorm(10)

temp <- log_normalize(log_x)

exp(temp[[2]]) - temp[[1]]

sum(exp(temp[[2]]))

log_x - temp[[3]] - temp[[2]]

}
