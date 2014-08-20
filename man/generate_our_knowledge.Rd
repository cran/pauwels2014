\name{generate_our_knowledge}
\alias{generate_our_knowledge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Initialize a knowledge list.
}
\description{
	This function defines a new knowledge list. Its main purpose is to give a value to all the global parameters that have to be set up before running the simulations.
}
\usage{
generate_our_knowledge(transform_params, global_parameters, datas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transform_params}{
	A function that transform parameter values. See \code{\link{transform_params}}. 
}
  \item{global_parameters}{
	A list of global parameters. See the details section for the default value.
}
  \item{datas}{
Potential data to be put in the \code{datas} slot.
}
}
\details{
There are four main slots in the knowledge list the first one is \code{transform_params} which is by default a log transformation function that allows to work in log space. See \code{\link{transform_params}}. The second slot, \code{global_parameters} contains various global parameters necessary for the simulation to proceed. The default value of those parameters is given in the next section. See also the vignette for details. The third slot is named \code{datas}. It is empty by default, but is intended to receive datasets, experiment function, posterior sample and expected risk estimates for each step of the simulation. The vignette provides more details about this.
}
\value{
A list that agregates the input. The default value of the global parameter is as follows:
\item{atol = 1e-6}{\code{deSolve} parameter (see \code{deSolve})}
\item{rtol = 1e-6}{\code{deSolve} parameter (see \code{deSolve})}
\item{tspan = seq(0,100,0.5)}{\code{deSolve} parameter, the time points at which all simulated kinetic quantities should be evaluated}
\item{max_step = 50}{\code{deSolve} parameter (see \code{deSolve})}
\item{max_it = 200}{optimization parameter, the maximum number of iteration in the for \code{\link{BFGS_special}} function}
\item{tol = 1e-3}{optimization parameter, a stoping criterion based on gradient norm tolerance for \code{\link{BFGS_special}} function}
\item{beta = 2}{optimization parameter, the division to perform on the step size when armijo's rule fails in \code{\link{BFGS_special}} function}
\item{c = 0.0001}{optimization parameter, armijo's rule second parameter for \code{\link{BFGS_special}} function}
\item{n_multi_mod_weight = 20}{sampling parameter, number of times the local search/sample operations are repeated in \code{\link{sample_function_multi_mod_weight}} function}
\item{max_log_like = - 700}{sampling parameter, a lower bound under which parameter values are systematically after the BFGS search in \code{\link{sample_function_multi_mod_weight}} function and in \code{\link{sample_function_single_mod}} function}
\item{centrality_ratio = 0.4}{sampling parameter, allows to only keep parameter values leading to reasonable datafits by filtering those which lead to kinetics trajectories that do not pass in the "middle" of the data in \code{\link{sample_function_multi_mod_weight}} function and in \code{\link{sample_function_single_mod}} function}
\item{sample_burn_in = 5000}{sampling parameter, size of the burn in sample to be ignored in the Metropolis Hasting algorithm algorithm in \code{\link{generate_sample}}}
\item{sample_to_keep1 = 10000}{sampling parameter, sample size to be further sampled in MH algorithm in \code{\link{generate_sample}}}
\item{sample_step = 1}{sampling parameter, the mean MH step length in \code{\link{generate_sample}} function}
\item{final_sample = 10000}{sampling parameter, the final sample size provided by sampling functions \code{\link{sample_function_multi_mod_weight}} and \code{\link{sample_function_single_mod}}. This is the size of the sample that is to be stored in the knowledge list.}
\item{final_sample_design = 100}{size of a subsample to be used for risk estimation in \code{\link{active_design}} strategy.}
\item{n_simu_weights = 100}{number of noise simulations required to estimate the weights in risk estimation in \code{\link{active_design}} strategy.}
\item{initial_conditions = c(g6 = 1, p6 = 1,p7 = 1,p8 = 1, v6_mrna = 0,v7_mrna = 0,v8_mrna = 0)}{simulation parameter, the initial condition default value}
\item{n_params = 9}{simulation parameter, number of free parameters}
\item{param_names =  c("p_degradation_rate", "r6_Kd", "r11_Kd", "pro6_strength", 
"pro7_strength", "pro9_strength", "rbs6_strength", "rbs7_strength", "rbs8_strength")}{simulation parameters, the names the free parameters}
\item{params = c(p_degradation_rate = 1, r6_Kd = 1, r11_Kd=1, pro6_strength = 1, pro7_strength = 1, pro9_strength = 1, rbs6_strength = 1, rbs7_strength = 1, rbs8_strength = 1)}{simulation parameters, an instance of free parameter named numeric vector}
\item{true_params = c(mrna6_degradation_rate =1, mrna7_degradation_rate =1, mrna8_degradation_rate =1, p_degradation_rate = 0.1,  r6_Kd = 2.6, r6_h = 4, r11_Kd=2, r11_h = 2, r12_Kd = 0.2, r12_h = 2, 
pro6_strength = 1,  pro7_strength = 0.8, pro9_strength = 3.77, rbs6_strength = 5,  rbs7_strength = 5, rbs8_strength = 5)}{simulation parameters, the true parameter to be estimated and used for simulations in the physical space}
\item{true_params_T = c(p_degradation_rate = 50, r6_Kd = 56.9162224661803, r11_Kd = 55.0171665943997, pro6_strength = 50, pro7_strength = 48.3848331165323, pro9_strength = 59.6056891700965, rbs6_strength = 61.649500072267, rbs7_strength = 61.649500072267, rbs8_strength = 61.649500072267)}{simulation parameters, the true parameter to be estimated and used for simulations in log space, such that \code{transform_params(true_params_T) = true_params}},
\item{dllname = "pauwels2014"}{simulation parameters, the name of the shared object which contains the function to be bassed to \code{deSolve} solver. See the vignette and \code{deSolve} package for more details. The default is \code{"pauwels2014"} as the source for the example of the article is provided in this package.}
}
\author{
Edouard Pauwels
}
\examples{
knobj <- generate_our_knowledge(transform_params)
knobj
}
\keyword{Model specific functions}

