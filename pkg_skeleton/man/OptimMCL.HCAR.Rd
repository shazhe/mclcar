\name{OptimMCL.HCAR}
\alias{OptimMCL.HCAR}

\title{Iterative procedure for maximising the Monte Carlo likelihood of the hierachical conditional auto-regressive 
models.}
\description{The function uses an iterative procedure of directly maximising the
Monte Carlo likelihood of a hierechical conditional auto-regressive model 
and the updating step size is limited by defining
an experimental region using the estimated Monte Carlo variance.}

\usage{
OptimMCL.HCAR <- function(data, psi0,  control = list())
}

\arguments{
\item{data}{A list or an environment contains the variables same as described in \code{\link{sim.HCAR}}.}
  \item{psi}{Starting value for the importance sampler parameter same as described in \code{\link{sim.HCAR}}.}
  \item{conrol}{a list of tuning parameters to control the algorithm. Details to be found at \code{\link{OptimMCL}}}
}

\value{
Same as in \code{\link{OptimMCL}}.
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{mcl.HCAR}}, \code{\link{sim.HCAR}}, \code{\link{summary.OptimMCL.HCAR}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CAR, HCAR}
\keyword{Monte Carlo likelihood}
