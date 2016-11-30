## The Scale can be fixed after some pilot runs for tuning
postZ <- function(data, family, psi, Z.start, mcmc.control=list(), plots = FALSE){
    N.data <- length(Z.start)
    mcmc.con <- list(method = "mala", N.Zy = 1e4, Scale = 1.65/(N.data^(2/6)), thin = 5, burns = 5e3,
                     scale.fixed = TRUE, c.ad = c(1, 0.7))
    ##namc <- names(mcmc.con)
    ##mcmc.con[(namc <- names(mcmc.control))] <- mcmc.control
    mcmc.con[names(mcmc.control)] <- mcmc.control


    switch(family,
           binom = postZ.binomial(data, Z.start, psi, mcmc.pars=mcmc.con, plots = plots),
           poisson = postZ.poisson(data, Z.start, psi, mcmc.pars=mcmc.con, plots = plots))
}


postZ.binomial <-function(data, Z.start, psi, mcmc.pars, plots = FALSE){
    method <- mcmc.pars$method
    switch(method,
           rwmh.iid = postZ.rwmh.binom(data, Z.start, psi, mcmc.pars, plots),
           rwmh.car = postZ.CAR.binom(data, Z.start, psi, mcmc.pars, plots),
           rwmh.Lap = postZ.Laplace.binom(data, Z.start, psi, mcmc.pars, plots),
           mala = postZ.mala.binom(data, Z.start, psi, mcmc.pars, plots))
}


postZ.poisson <-function(data, Z.start, psi, mcmc.pars, plots = FALSE){
    method <- mcmc.pars$method
    switch(method,
           rwmh.iid = postZ.rwmh.poisson(data, Z.start, psi, mcmc.pars, plots),
           rwmh.car = postZ.CAR.poisson(data, Z.start, psi, mcmc.pars, plots),
           rwmh.Lap = postZ.Laplace.poisson(data, Z.start, psi, mcmc.pars, plots),
           mala = postZ.mala.poisson(data, Z.start, psi, mcmc.pars, plots))
}

################################################################
####     Simulate data from the joint distribution of       ####
####     the observed data Y and latent CAR variable Z      ####
####     for binomial or Poisson data glm model             ####
################################################################
CAR.JointGLM <- function(method = c("binom", "poisson"), N, data, psi){
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    n.trial <- data$n.trial
    W <- data$W
    covx <- data$covX
    XB <- covx %*% beta

    n <-nrow(W)
    YZ <- rep(list(0),N)

    I <- diag.spam(1, n)
    Q <- as.spam((I - rho*W)/sigma)
    L <- chol(Q,pivot="MMD",  memory = list(nnzcolindices = 6.25 * n))
    if(method == "binom"){
        for (i in 1:N){
            zz<- solve.spam(L, rnorm(n))
            eta <- XB + zz
            p <- as.vector(exp(eta)/(1+exp(eta)))
            yy <- rbinom(n, n.trial, p)
            YZ[[i]] <- cbind(yy,zz)
        }}else{
            for (i in 1:N){
            zz<- solve.spam(L, rnorm(n))
            eta <- XB + zz
            lambda <- exp(eta)
            yy <- rpois(n, lambda)
            YZ[[i]] <- cbind(yy,zz)
        }}
    return(YZ)
}
