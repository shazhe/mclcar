\name{postZ}
\alias{postZ}
\title{
Sampling the CAR latent variables given the Binomial or Poisson
observations in glm models.
}
\description{
The function uses several different MCMC algorithms to sample the CAR
latent variables Z from the distribution P(Z|y) given the Binomial or
Poisson observations in glm models using the importance sampler parameter
values psi.
}
\usage{
postZ(data, family, psi, Z.start, mcmc.control = list(), plots = FALSE)
}
\arguments{
  \item{data}{ the data set to be estimated and the object has same as
	the output of \code{\link{CAR.simGLM}}
}
  \item{family}{ a character take values in "binom" and "poisson"
	for Binomial or Poisson glm with CAR latent variable).
}
  \item{psi}{
Parameters for the importance distribution.
}
  \item{Z.start}{
Initial value of the CAR variables for the Markov Chain.
}
  \item{mcmc.control}{
a list that controls the MCMC algorithm, see more in Details.
}
  \item{plots}{
if TRUE, the diagnostic plots are shown on running the function.
}
}
\details{The function uses two MCMC algorithms, the Metropolis Hastings
  algorithm with random walk proposal (rwmh) and the Metropolis adjusted
  Langevin algorithm (mala), in sampling the CAR latent variables.
  The scale parameter for the proposal distribution can be tuned to
  achieve the optimal acceptance rate by using an adaptive algorithm in
  a pilot run. The argument \code{mcmc.control} contains the following
  objects
  \itemize{
	\item{N.Zy,}{ the number of iterations of the Markov Chain}
	\item{Scale,}{ the scale for the proposal distribution}
	\item{thin,}{ the thinning step}
	\item{burns,}{ the number of burn-in samples to be discarded}
	\item{method,}{ the name of the MCMC algorithm, either "rwmh" or "mala"}
	\item{scale.fixed,}{ if true, the scale parameter is fixed;
	  otherwise an adaptive MCMC algorithm is run for tuning the scale
	  parameter to achieve the optimal acceptance rate by using an
	  stochastic approximation algorithm}
	\item{c.ad,}{ an array contains the two values for the tuning
	  parameter for the stochastic approximation algorithm in the
	  adaptive MCMC algorithm; default is set to be c(1, 0.7)}
	}
}
\value{
  The function return a list of the following objects:
  \itemize{
	\item{Z} {a matrix of samples from from P(Z|y)}
	\item{AC.r}{an array of the acceptance rate updated at each
	  iteration}
	\item{Scales}{if in the list of mcmc.control scale.fixed = FALSE, an
	array of the adaptive scales for the MCMC proposal}
  }
}

\author{
Zhe Sha \email{sha@stats.ox.ac.uk}
}

\seealso{
\code{\link{mcl.prep.glm}}, \code{\link{OptimMCL}}, \code{\link{rsmMCL}}
}
\examples{
set.seed(33)
n.torus <- 10
nb <- 30
rho <- 0.2
sigma <- 1.5
beta <- c(1, 1)
pars.true <- c(rho, sigma, beta)
X0 <- cbind(rep(1, n.torus^2), sample(log(1:n.torus^2)/5))
mydata2 <- CAR.simGLM(method = "binom", n = c(n.torus, n.torus), pars = pars.true,
                      Xs = as.matrix(X0), n.trial = nb)

## use a glm to find initial values for the importance sampler
library(spdep)
data.glm <- data.frame(y=mydata2$y, mydata2$covX[,-1])
fit.glm <- glm(cbind(y, nb-mydata2$y) ~ .,data = data.glm, family=binomial)
## estimate sigma and rho, transform the binomial to Gaussian by logit
logitp <- log((mydata2$y+0.5)/(mydata2$n.trial - mydata2$y + 0.5))
data.splm <- data.frame(y=logitp, mydata2$covX[,-1])
listW <- mat2listw(mydata2$W)
fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
pars1 <- c(fit.splm$lambda, fit.splm$fit$s2, coef(fit.glm))

## Sample form importance distribution with psi = pars1
mc.control <- list(N.Zy = 1e3, Scale = 1.65/(n.torus^(2/6)), thin = 5,
                   burns = 5e2, method = "mala", scale.fixed = TRUE)
## Binomial
Z.S0 <- CAR.simWmat(pars1[1], 1/pars1[2], mydata2$W)
simZy <- postZ(data = mydata2, Z.start = Z.S0, psi = pars1,
               family = "binom", mcmc.control = mc.control, plots =
TRUE)

}

\keyword{MCMC}
\keyword{Spatial}
