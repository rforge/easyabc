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
  \item{method}{The choosen method for performing the sequential ABC}
  \item{model}{The function to invoke for running the model}
  \item{prior_matrix}{A matrix containing the range of the prior distribution of the parameters}
  \item{nb_simul}{The number of simulations}
  \item{tolerance_tab}{An array containing the tolerance levels}
  \item{summary_stat_target}{An array containing the summary statistics of targeted (observed) data}
  \item{\dots}{The additional arguments depending of the choosen method}
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
%% ~put references to the literature/web site here ~
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
ABC_sequential(method="Beaumont",model=trait_model,prior_matrix=priormatrix,
nb_simul=10,summary_stat_target=sum_stat_obs,tolerance_tab=c(8,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{abc}

