\name{transform_params}
\alias{transform_params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User defined parameter transformation function.
}
\description{
Allows to go from log space to physical space for kinetic parameters. Known parameter values are also added here. 
}
\usage{
transform_params(pars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
A named numeric parameter vector
}
}
\details{
The number and order of resulting parameters should be the same as the one expected by the \code{\link{simulate_experiment_no_transform}} function.
}
\value{
A named numeric parameter vector
}
\author{
Edouard Pauwels
}

\seealso{
\code{\link{simulate_experiment}}, \code{\link{simulate_experiment_no_transform}}, \code{\link{estimate_risk_out_all}}
}
\examples{
## Generate the knowledge object with correct parameter value
knobj <- generate_our_knowledge(transform_params)

knobj$global_parameters$true_params_T
transform_params(knobj$global_parameters$true_params_T)
}
\keyword{Model specific functions}
