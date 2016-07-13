###############################################################################################
####                              Main Evaluation Functions                                ####
###############################################################################################
#### preparing MCL samples.
mcl.prep.glm <- function(data, family, psi, Z.start = NULL,
                         mcmc.control = list(), pilot.run = TRUE, pilot.plot = FALSE, plot.diag = FALSE){
    if(is.null(Z.start)){
        Z.start <- CAR.simWmat(psi[1], 1/psi[2], data$W)
    }

    N.data <- length(Z.start)
    mcmc.pars <- list(method = "mala", N.Zy = 1e4, Scale = 1.65/(N.data^(2/6)), thin = 5, burns = 5e3,
                     scale.fixed = TRUE, c.ad = c(1, 0.7))

    mcmc.pars[names(mcmc.control)] <- mcmc.control

    if(pilot.run == TRUE){
        ad.pars<- mcmc.pars
        ad.pars$N.Zy <- min(5e3, mcmc.pars$N.Zy/10)
        ad.pars$burns <- 0
        ad.pars$scale.fixed <- FALSE
        Zy <- postZ(data = data, family, psi, Z.start, mcmc.control = ad.pars, plots = pilot.plot)

        Scale.fix <- mean(Zy[[3]][(ad.pars$N.Zy/5+1):(ad.pars$N.Zy)])
        mcmc.pars$Scale <- Scale.fix
        Z.start <- Zy[[1]][ad.pars$N.Zy/ad.pars$thin, ]
        Zy <- postZ(data = data, family, psi, Z.start, mcmc.control = mcmc.pars, plots = plot.diag)
    }
    else{
         Zy <- postZ(data = data, family, psi, Z.start, mcmc.control = mcmc.pars, plots = plot.diag)
   }

    ## Calculate the fixed parameter in the approximate
    Zys <- Zy[[1]]
    logH.Zy.psi <- apply(Zys, 1, logH.psi, pars = psi, data = data, family = family)
    data$Zy <- Zys
    data$lHZy.psi <- logH.Zy.psi
    data$mcmc.pars <- mcmc.pars
    return(data)
}

#### Monte Carlo likelihood
mcl.glm <- function(pars, mcdata, family, Evar = FALSE){
    switch(family,
           binom = mcl.bin(pars, mcdata, Evar),
           poisson = mcl.pois(pars, mcdata, Evar))
}

#### profiling MC-MLE of beta
get.beta.glm <- function(beta0, rho.sig , mcdata, family){
    switch(family,
           binom = get.beta.bin(beta0, rho.sig, mcdata),
           poisson = get.beta.pois(beta0, rho.sig, mcdata))
}

#### Monte Carlo profile likelihood
mcl.profile.glm <- function(pars, beta0, mcdata, family, Evar = FALSE){
    switch(family,
           binom = mcl.bin.profile(pars, beta0, mcdata, Evar),
           poisson = mcl.pois.profile(pars, beta0, mcdata, Evar))
}

#### Monte Carlo Variance
vmle.glm <- function(MLE, mcdata, family){
    switch(family,
           binom = vmle.bin(MLE, mcdata),
           poisson = vmle.pois(MLE, mcdata))
}

#### Unormalised density, gradient and Hessian
logH.psi <- function(pars, data, mcdata, family){
    switch(family,
           binom = logH.psi.bin(pars, data, mcdata),
           poisson = logH.psi.pois(pars, data, mcdata))
}


mcl.grad.glm<- function(pars, mcdata, family, Evar = FALSE){
    switch(family,
           binom = mcl.grad.bin(pars, mcdata, Evar),
           poisson = mcl.grad.pois(pars, mcdata, Evar))
}


mcl.Hessian.glm <- function(pars, mcdata, family){
    switch(family,
           binom = mcl.Hessian.bin(pars, mcdata),
           poisson = mcl.Hessian.pois(pars, mcdata))
}
