\name{sample_function}
\alias{sample_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates posterior samples
}
\description{
This generates a posterior sample based on the data and global parameters defined in its input.
}
\usage{
sample_function(knobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knobj}{
A knowledge list. See \code{\link{knobjs}}.
}
}
\details{
This function is the same as, \code{\link{sample_function_multi_mod_weight}}. This is the sampling function used by design strategies. It can be changed to modify the sampling strategy (e.g. \code{\link{sample_function_single_mod}}).
}
\value{
A matrix which row represent a parameter named numeric vector in log space. Its distribution is supposed to be closed to the posterior distribution defined by knobj.
}

\author{
Edouard Pauwels
}

\seealso{
\code{\link{sample_function_multi_mod_weight}}, \code{\link{sample_function_single_mod}}
}

\keyword{Sampling functions}

