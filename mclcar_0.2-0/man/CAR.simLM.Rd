\name{CAR.simLM}
\alias{CAR.simLM}
\alias{CAR.simGLM}
\alias{CAR.simTorus}
\alias{CAR.simWmat}

\title{Simulate samples from a CAR model.
}
\description{This help page documents several functions for simulate
  samples from different CAR models.
}

\usage{
CAR.simTorus(n1, n2, rho, prec)
CAR.simWmat(rho, prec, W)
CAR.simLM(pars, data)
CAR.simGLM(method = c("binom", "poisson"), W, n = NULL, pars, Xs = NULL,
n.trial = 1)

}

\arguments{
  \item{n1, n2, n}{n = c(n1, n2), the size of the torus}
  \item{rho}{the spatial coefficient in the CAR precision matrix (inverse of the covariance matrix)}
  \item{prec}{the precision in the CAR precision matrix}
  \item{W}{the spatial weight matrix}
  \item{pars}{parameter values of the direct CAR model}
  \item{data}{a \code{data} object same as described in \code{\link{loglik.dCAR}}}
  \item{method}{ a character equal to either "binom" or "poisson"
	that indicate the distribution of the sample to be simulated.}
  \item{Xs}{the covariates to be used in the glm}
  \item{n.trial}{the number of trials in each unit of the binomial samples}
}

\value{
  \item{\code{CAR.simTorus}}{returns a list containing the spatial
	weight matrix \code{W} and the simulated CAR samples \code{X};}
  \item{\code{CAR.simWmat},\code{CAR.simLM} }{returns an array of the
  simulated CAR samples;}
\item{\code{CAR.simGLM}}{returns a list containing:
  \describe{
	\item{rho, sigma, beta,}{the parameter values used to generate the
	  samples}
	\item{y,}{the simulated data}
	\item{covX,}{the covariates}
	\item{W,}{the spatial weight matrix}
	\item{Z.true,}{the simulated CAR variables}
	\item{eta,}{the simulated glm linear response}
	\item{Emean,}{the simulated glm mean}
	\item{n.trial,}{the number of trials in each unit for binomial samples}
	}}
}

\author{Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{mcl.glm}}, \code{\link{mcl.dCAR}},
\code{\link{mcl.prep.dCAR}}, \code{\link{mcl.prep.glm}}
}
\examples{
## Simulate CAR data on a torus
set.seed(33)
n.torus <- 10
rho <- 0.2
sigma <- 1.5
prec <- 1/sigma
beta <- c(1, 1)
XX <- cbind(rep(1, n.torus^2), sample(log(1:n.torus^2)/5))
mydata1 <- CAR.simTorus(n1 = n.torus, n2 = n.torus, rho = rho, prec = prec)

## Simulate CAR data for a given spatial weight matrix
Wmat <- mydata1$W
mydata2 <- CAR.simWmat(rho = rho, prec = prec, W = Wmat)

## Simualte data from a linear model with CAR error
y <- XX \%*\% beta + mydata1$X
mydata1$data.vec <- data.frame(y=y, XX[,-1])
mydata3 <- CAR.simLM(pars = c(0.1, 1, 2, 0.5), data = mydata1)

## Simulate Binomial data with CAR latent variables
mydata4 <- CAR.simGLM(method="binom", n=c(10,10), pars = c(rho, sigma,
                      beta), Xs=XX, n.trial = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{CAR}
\concept{Monte Carlo likelihood}
