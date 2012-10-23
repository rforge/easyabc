\name{ABC_mcmc_cluster}
\alias{ABC_mcmc_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ABC_mcmc_cluster(method, model, prior_matrix, n_obs, n_between_sampling, summary_stat_target, n_cluster = 1, ...)
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
  \item{n_obs}{
%%     ~~Describe \code{n_obs} here~~
}
  \item{n_between_sampling}{
%%     ~~Describe \code{n_between_sampling} here~~
}
  \item{summary_stat_target}{
%%     ~~Describe \code{summary_stat_target} here~~
}
  \item{n_cluster}{
%%     ~~Describe \code{n_cluster} here~~
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
\author{
%%  ~~who you are~~
}
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
function (method, model, prior_matrix, n_obs, n_between_sampling, 
    summary_stat_target, n_cluster = 1, ...) 
{
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior_matrix)) 
        stop("'prior_matrix' is missing")
    if (missing(n_obs)) 
        stop("'n_obs' is missing")
    if (missing(n_between_sampling)) 
        stop("'n_between_sampling' is missing")
    if (missing(summary_stat_target)) 
        stop("'summary_stat_target' is missing")
    if (!any(method == c("Marjoram_original", "Marjoram", "Wegmann"))) {
        stop("Method must be Marjoram_original, Marjoram or wegmann")
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
    if (!is.vector(n_cluster)) 
        stop("'n_cluster' has to be a number.")
    if (length(n_cluster) > 1) 
        stop("'n_cluster' has to be a number.")
    if (n_cluster < 1) 
        stop("'n_cluster' has to be a positive number.")
    n_cluster = floor(n_cluster)
    options(scipen = 50)
    switch(EXPR = method, Marjoram = .ABC_MCMC2_cluster(model, 
        prior_matrix, n_obs, n_between_sampling, summary_stat_target, 
        n_cluster, ...), Wegmann = .ABC_MCMC3_cluster(model, 
        prior_matrix, n_obs, summary_stat_target, n_cluster, 
        ...))
    options(scipen = 0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
