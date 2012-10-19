\name{ABC_rejection}
\alias{ABC_rejection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{brute-force ABC (Pritchard et al. 1999)}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ABC_rejection(model, prior_matrix, nb_simul, use_seed = TRUE, seed_count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{prior_matrix}{The matrix of the priors coefficients}
  \item{nb_simul}{number of time to evaluate the model}
  \item{use_seed}{
%%     ~~Describe \code{use_seed} here~~
}
  \item{seed_count}{
%%     ~~Describe \code{seed_count} here~~
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
[Pritchard et al. 1999] J.K. Pritchard, M.T. Seielstad, A. Perez-Lezaun, M.W. Feldman (1999) Population growth of human Y chromosomes: a study of Y chromosome microsatellites., Mol. Biol. Evol., 16, 1791–1798 
}
\author{Franck Jabot, Thierry Faure, Nicolas Dumoulin}\note{
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
function (model, prior_matrix, nb_simul, use_seed = TRUE, seed_count = 0) 
{
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior_matrix)) 
        stop("'prior_matrix' is missing")
    if (missing(nb_simul)) 
        stop("'nb_simul' is missing")
    if (!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) 
        stop("'prior_matrix' has to be a matrix or data.frame.")
    if (is.data.frame(prior_matrix)) 
        prior_matrix <- as.matrix(prior_matrix)
    if (dim(prior_matrix)[2] != 2) 
        stop("'prior_matrix' must have two columns.")
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1.")
    nb_simul = floor(nb_simul)
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean.")
    if (!is.vector(seed_count)) 
        stop("'seed_count' has to be a number.")
    if (length(seed_count) > 1) 
        stop("'seed_count' has to be a number.")
    if (seed_count < 0) 
        stop("'seed_count' has to be a positive number.")
    seed_count = floor(seed_count)
    options(scipen = 50)
    tab_simul_summarystat = NULL
    tab_param = NULL
    start = Sys.time()
    pb <- .progressBar(width = 50)
    duration = 0
    for (i in 1:nb_simul) {
        l = dim(prior_matrix)[1]
        param = array(0, l)
        for (j in 1:l) {
            param[j] = runif(1, min = prior_matrix[j, 1], max = prior_matrix[j, 
                2])
        }
        if (use_seed) {
            param = c((seed_count + i), param)
        }
        simul_summarystat = model(param)
        tab_simul_summarystat = rbind(tab_simul_summarystat, 
            simul_summarystat)
        if (use_seed) {
            tab_param = rbind(tab_param, param[2:(l + 1)])
        }
        else {
            tab_param = rbind(tab_param, param)
        }
        duration = difftime(Sys.time(), start, unit = "secs")
        text = ""
        if (i == nb_simul) {
            text = paste("Completed  in", format(.POSIXct(duration, 
                tz = "GMT"), "\%H:\%M:\%S"), "                                              ")
        }
        else {
            text = paste("Time elapsed:", format(.POSIXct(duration, 
                tz = "GMT"), "\%H:\%M:\%S"), "Estimated time remaining:", 
                format(.POSIXct(duration/i * (nb_simul - i), 
                  tz = "GMT"), "\%H:\%M:\%S"))
        }
        .updateProgressBar(pb, i/nb_simul, text)
    }
    close(pb)
    options(scipen = 0)
    sd_simul = sapply(as.data.frame(tab_simul_summarystat), sd)
    list(param = tab_param, stats = tab_simul_summarystat, weights = array(1/nb_simul, 
        nb_simul), stats_normalization = sd_simul, nsim = nb_simul, 
        computime = as.numeric(difftime(Sys.time(), start, unit = "secs")))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
