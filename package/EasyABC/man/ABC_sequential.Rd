\name{ABC_sequential}
\alias{ABC_sequential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sequential ABC methods}
\description{Sequential ABC methods 
It contains algorithms from Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)}
\usage{
ABC_sequential(functionName, model, prior_matrix, nb_simul, tolerance_tab, summary_stat_target, use_seed = TRUE, seed_count = 0, inside_prior = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{functionName}{
%%     ~~Describe \code{functionName} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{prior_matrix}{
%%     ~~Describe \code{prior_matrix} here~~
}
  \item{nb_simul}{
%%     ~~Describe \code{nb_simul} here~~
}
  \item{tolerance_tab}{
%%     ~~Describe \code{tolerance_tab} here~~
}
  \item{summary_stat_target}{
%%     ~~Describe \code{summary_stat_target} here~~
}
  \item{use_seed}{
%%     ~~Describe \code{use_seed} here~~
}
  \item{seed_count}{
%%     ~~Describe \code{seed_count} here~~
}
  \item{inside_prior}{
%%     ~~Describe \code{inside_prior} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(functionName,model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=TRUE,seed_count=0,inside_prior=TRUE){
	## fonction ? faire qui regroupe les algos sequentiels (PMC, drovandi, delmoral, maxime)
  switch(EXPR = functionName, 
         PMC = .ABC_PMC(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
         PMC2 = .ABC_PMC2(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
         Drovandi = .ABC_Drovandi(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha,c,first_tolerance_level_auto,use_seed,seed_count),
         Delmoral = .ABC_Delmoral(model,prior_matrix,nb_simul,alpha,M,nb_threshold,tolerance_target,summary_stat_target,use_seed,seed_count),
         Lenormand = .ABC_Lenormand(model,prior_matrix,nb_simul,summary_stat_target,alpha,p_acc_min,use_seed,seed_count,inside_prior))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
