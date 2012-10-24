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

prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
ABC_rejection(model=trait_model, prior_matrix=prior_matrix, nb_simul=10, use_seed=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ abc }
