\name{ABC_sequential}
\alias{ABC_sequential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sequential ABC methods}
\description{Sequential ABC methods 
It contains algorithms from Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)}
\usage{
ABC_sequential(method, model, prior_matrix, nb_simul, tolerance_tab, summary_stat_target, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{The choosen method for performing the sequential ABC. It can take the following character values:"Beaumont", "Drovandi", "Delmoral", "Lenormand"}
  \item{model}{The function to invoke for running the model}
  \item{prior_matrix}{A matrix containing the range of the prior distribution of the parameters}
  \item{nb_simul}{The number of simulations to perform}
  \item{summary_stat_target}{An array containing the summary statistics of targeted (observed) data}
  \item{use_seed}{use of the seed value. Default value TRUE}
  \item{seed_count}{initial value for seed number. Default value 0}
  \item{verbose}{for debugging purpose. default value FALSE}
  \item{inside_prior}{Be carefull this parameter is only used by Beaumont and Lenormand algorithms.Default value TRUE for Beaumont and FALSE for Lenormand.}
  \item{tolerance_tab}{An array containing the tolerance levels. Default value -1. Be carefull this parameter is only used by Beaumont and Drovandi algorithms.}
  \item{alpha}{Default value 0.5. Be carefull this parameter is only used by Drovandi, Delmoral and Lenormand algorithms.}
  \item{c}{Default value 0.01. Be carefull this parameter is only used by Drovandi algorithm.}
  \item{first_tolerance_level_auto}{Default value TRUE. Be carefull this parameter is only used by Drovandi algorithm.}
  \item{M}{Be carefull this parameter is only used by Delmoral algorithm.}
  \item{nb_threshold}{Be carefull this parameter is only used by Delmoral algorithm.}
  \item{tolerance_target}{Default value 1. Be carefull this parameter is only used by Delmoral algorithm.}
  \item{p_acc_min}{Default value 0.05. Be carefull this parameter is only used by Lenormand algorithm.}
}
\details{
If a blank string is given as \code{method}, the list of available methods will be printed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
    [Beaumont et al., 2009] Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika,96(4):983-990.
    [Drovandi & Pettitt 2011] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. Biometrics, 67(1):225-233.
    [Del Moral et al. 2012] Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation, Statistics and Computing., 22(5):1009-1020.
   [Lenormand et al. 2012] Lenormand, M., Jabot, F., Deffuant G. (2012). Adaptive approximate Bayesian computation for complex models, submitted to Comput. Stat. )
}
\author{Franck Jabot, Thierry Faure, Nicolas Dumoulin}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
priormatrix=c(500,3,-2.3,1,-25,-0.7,500,5,1.6,1,125,3.2)
dim(priormatrix)<-c(6,2)
sum_stat_obs=c(100,2.5,20,30000)
 ## to perform the Beaumont Method
ABC_sequential(method="Beaumont",model=trait_model,prior_matrix=priormatrix,
nb_simul=10,summary_stat_target=sum_stat_obs,tolerance_tab=c(8,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{abc}

