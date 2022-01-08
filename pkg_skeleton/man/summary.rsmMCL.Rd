\name{summary.rsmMCL}
\alias{summary.rsmMCL}
\title{Summary the output from the response surface method of maximising
  the Monte Carlo likelihood
}
\description{
This function summarizes the output of the output from the
  function \code{\link{rsmMCL}}.
}
\usage{
\method{summary}{rsmMCL}(object, family, trace.all = TRUE, mc.covar =
TRUE, \dots)
}

\arguments{
  \item{object}{ an rsmMCL object returned by \code{\link{rsmMCL}}.
}
  \item{family}{
a character takes value in "gauss", "binom" and "poisson".
}
  \item{trace.all}{
an logic value tells whether the input object given
  by \code{\link{rsmMCL}} contains results from all iterations of not
}
  \item{mc.covar}{
if TRUE, the estimated covariance matrix of the MC-MLE is returned
}
\item{\dots}{
arguments passed to or from other methods.}
}

\value{
   A list containing the following objects:
  \describe{
	\item{MC.mle,}{ the final MC-MLE}
	\item{N.iter,}{ the total number of iterations}
	\item{total.time,}{ the total time elapsed}
	\item{convergence,}{ if TRUE the procedure converges}
	\item{hessian,}{ the Hessian at the MC-MLE if given; the default is NULL}
	\item{mc.covar}{ the estimated covariance matrix of the MC-MLE if
  given; the default is NULL}
	\item{mc.samples}{the Monte Carlo samples size used in the initial
  stage and after the first convergence.}
  }
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{rsmMCL}}
}
\examples{
## See examples for rsmMCL
}

\concept{CAR}
\concept{Monte Carlo likelihood}
\concept{Spatial}
\concept{Response surface design}