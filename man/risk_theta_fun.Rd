\name{risk_theta_fun}
\alias{risk_theta_fun}

\title{
Risk function
}
\description{
Compute the risk between two different parameter values.
}
\usage{
risk_theta_fun(theta1, theta2, n_params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta1}{
A first parameter array.
}
  \item{theta2}{
A second parameter array.
}
  \item{n_params}{
Number of parameters to be estimated in the two arrays.
}
}
\details{
The default is averaged squared difference between log values of parameters 
}
\value{
A numerical value.
}
\author{
Edouard Pauwels
}
\examples{
theta1 <- sample(1:10000, size = 100)
theta2 <- sample(1:10000, size = 100)
n <- 50

risk_theta_fun(theta1, theta2, n)
}
\keyword{Model specific functions}
