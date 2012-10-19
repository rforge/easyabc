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
  \item{method}{
%%     ~~Describe \code{method} here~~
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
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (method, model, prior_matrix, nb_simul, summary_stat_target, 
    ...) 
{
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior_matrix)) 
        stop("'prior_matrix' is missing")
    if (missing(nb_simul)) 
        stop("'nb_simul' is missing")
    if (missing(summary_stat_target)) 
        stop("'summary_stat_target' is missing")
    if (!any(method == c("Beaumont", "Drovandi", "Delmoral", 
        "Lenormand"))) {
        stop("Method must be Beaumont, Drovandi, Delmoral or Lenormand")
    }
    if (!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) 
        stop("'prior_matrix' has to be a matrix or data.frame.")
    if (is.data.frame(prior_matrix)) 
        prior_matrix <- as.matrix(prior_matrix)
    if (dim(prior_matrix)[2] != 2) 
        stop("'prior_matrix' must have two columns.")
    if (!is.vector(nb_simul)) 
        stop("'nb_simul' has to be a number.")
    if (length(nb_simul) > 1) 
        stop("'nb_simul' has to be a number.")
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1.")
    nb_simul = floor(nb_simul)
    if (!is.vector(summary_stat_target)) 
        stop("'summary_stat_target' has to be a vector.")
    options(scipen = 50)
    switch(EXPR = method, Beaumont = .ABC_PMC(model, prior_matrix, 
        nb_simul, summary_stat_target, ...), Drovandi = .ABC_Drovandi(model, 
        prior_matrix, nb_simul, summary_stat_target, ...), Delmoral = .ABC_Delmoral(model, 
        prior_matrix, nb_simul, summary_stat_target, ...), Lenormand = .ABC_Lenormand(model, 
        prior_matrix, nb_simul, summary_stat_target, ...))
    options(scipen = 0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
