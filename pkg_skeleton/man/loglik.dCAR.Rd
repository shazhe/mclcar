\name{loglik.dCAR}
\alias{loglik.dCAR}
\alias{ploglik.dCAR}
\alias{get.beta.lm}
\alias{sigmabeta}
\alias{mple.dCAR}

\title{Likelihood computing and parameter estimation for a direct CAR model.
}
\description{This help page documents a few functions for exact evaluation of the
  log-likelihood and profile log-likelihood for a direct CAR model, the
  maximum pseudo-likelihood estimator, and least square estimators for
  beta and sigma given the spatial coefficient rho in the CAR covariance
  matrix. The exact evaluation computes the log determinant in the log-likelihood
  function with eigen-values of the spatial weight matrix.
}

\usage{
loglik.dCAR(pars, data)
ploglik.dCAR(rho, data)

get.beta.lm(rho, data)
sigmabeta(rho, data)

mple.dCAR(data, tol = 1e-06, rho0=0)
}

\arguments{
  \item{pars}{the parameter value to be evaluated}
   \item{rho}{the value of the spatial coefficient in the CAR covariance
  to be evaluated}
\item{data}{a list containing the following objects:
\describe{
  \item{\eqn{W}}{the spatial weight matrix,}
  \item{lambda}{the eigenvalues of \eqn{W} (if given),}
  \item{data.vec}{a data frame of the response variable and covariates.}
}
}
  \item{tol}{tolerance for the relative difference for two consequtive
	iterations in finding the MPLE; default set to be = 1e-06}

  \item{rho0}{starting value for iteratively finding the MPLE; default
  value set to be 0.}
}

\details{
The eigen-values can be supplied in the \code{data} list if the
likelihood for the same data is going to be evaluated many times; if not
supplied the function use the function \code{\link{eigen}} to find the
eigen-values of the weight matrix.
}
\value{\code{loglik.dCAR} and \code{ploglik.dCAR} return a numeric value
  of the log-likelihood evaluation. \code{get.beta.lm}, \code{sigmabeta}
  and \code{mple} return a numeric array of the estimates.
}

\author{
Zhe Sha \email{sha@stats.ox.ac.uk}
}

\seealso{
\code{\link{CAR.simLM}}, \code{\link{get.beta.glm}}, \code{\link{mcl.dCAR}}
}
\examples{
## Simulate data from a torus
set.seed(30)
n.torus <- 20
rho <- 0.15
sigma <- 1.5
beta <- c(1, 2)
XX <- cbind(rep(1, n.torus^2), log(1:n.torus^2))
mydata <- CAR.simTorus(n.torus, n.torus, rho, 1/sigma)
y <- XX \%*\% beta + mydata$X
mydata$data.vec <- data.frame(y=y, XX[,-1])

## evaluate the log-likelihood without lambda
loglik.dCAR(pars = c(0.1, 1, 0.9, 2.1), data = mydata)

## evaluate the log-likelihood with lamda
lambda <- eigen(mydata$W, symmetric = TRUE, only.values=TRUE)$values
mydata$lambda <- lambda
loglik.dCAR(pars = c(0.1, 1, 0.9, 2.1), data = mydata)

## evaluate the profile log-likelihood of rho
ploglik.dCAR(rho = 0.1, data = mydata)

## given rho = 0.1, find the least square estimates for beta and sigma
get.beta.lm(rho = 0.1, data = mydata)
sigmabeta(rho = 0.1, data = mydata)

## find the maximum pseudo-likelihood estimates
mple.dCAR(data = mydata)
}
\keyword{CAR}
\keyword{likelihood}
