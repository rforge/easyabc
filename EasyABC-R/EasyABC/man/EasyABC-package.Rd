\name{EasyABC-package}
\alias{EasyABC-package}
\alias{EasyABC}
\docType{package}
\title{A package for ABC computation}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: EasyABC
Type: Package
Version: 1.0
Date: 2012-10-15
License: GPL-3
LazyLoad: yes
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{Franck Jabot, Thierry Faure, Nicolas Dumoulin}
\references{
in preparation
}
\keyword{ package, abc }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix
}