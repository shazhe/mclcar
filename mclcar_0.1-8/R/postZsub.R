###############################################################################################
####                               Sub Simulation Functions                                ####
###############################################################################################
#### Simulation from posterior from glm binomial with latent CAR
#### using different proposals
########################################################################################
postZ.rwmh.binom <-  function(data, Z.start, psi, mcmc.pars, plots = FALSE){
    #for(i in 1:length(data)){
    # tempobj=data[[i]]
    # eval(parse(text=paste(names(data)[[i]],"= tempobj")))
 #}
  #   for(i in 1:length(mcmc.pars)){
   #  tempobj=mcmc.pars[[i]]
    # eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
 #}
	W <- data$W
	covX <- data$covX
	y <- data$y
	n.trial <- data$n.trial
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    t0 <- sum(y * Z0) - sum(n.trial * log(1 + exp(xb + Z0))) - Z0 %*% Q %*% Z0 /2

    acc <- 0
    accept <- rep(0,N.Zy)

    if (scale.fixed){ ## fixed Scale
        for(i in 1:N.Zy) {
            Zb <- Z0 + rnorm(nn, sd = Scale)
            t1 <- sum(Zb * y) - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <-acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept))
    }
    else { ## adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <- Z0 + rnorm(nn, sd = Scale)
            t1 <- sum(Zb * y) - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <-acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept, Scale.vec))
    }
}


######################################################################################################
## Use RWHM with CAR proposal to sample from Z|y
postZ.CAR.binom <-  function(data, Z.start, psi, mcmc.pars, plots = FALSE){
     #for(i in 1:length(data)){
     #tempobj=data[[i]]
     #eval(parse(text=paste(names(data)[[i]],"= tempobj")))
 #}
     #for(i in 1:length(mcmc.pars)){
     #tempobj=mcmc.pars[[i]]
     #eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
# }
 	W <- data$W
	covX <- data$covX
	y <- data$y
	n.trial <- data$n.trial
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    L <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * nn)) ## upper triagle
    t0 <- y %*% Z0 - sum(n.trial * log(1 + exp(xb + Z0))) - Z0 %*% Q %*% Z0 / 2

    acc <- 0
    accept <- rep(0,N.Zy)

    ##if (scale.fix){
      if (scale.fixed){
      for(i in 1:N.Zy) { ## Fixed scale
            Zb <-Z0 + backsolve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc + 1
            }
            accept[i] <-acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept))
    }
    else{ ## Adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <-Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc + 1
            }
            accept[i] <-acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept, Scale.vec))
    }
}

#######################################################################################
## Use RWHM with Laplace approximation as proposal to sample from Z|y
postZ.Laplace.binom  <-  function(data, Z.start, psi, mcmc.pars, plots = FALSE){
    #for(i in 1:length(data)){
    # tempobj=data[[i]]
    # eval(parse(text=paste(names(data)[[i]],"= tempobj")))
 #}
  #   for(i in 1:length(mcmc.pars)){
  #   tempobj=mcmc.pars[[i]]
  #   eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
 #}
  	W <- data$W
	covX <- data$covX
	y <- data$y
	n.trial <- data$n.trial
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma

    ## Laplace approximation to the posterior target
    logpost.Lap <- function(z){
        I <- as.matrix(diag.spam(1, nrow(W)))
        t1 <- sum((z+xb) %*% y)
        t2 <- -sum(log(1+exp(z+xb)))
        t3 <- -0.5 * z %*% (I - rho*W) %*% z
        sum(t1,t2,t3)
    }
    dlogpost.Lap <- function(z){
        I <- as.matrix(diag.spam(1, nrow(W)))
        dz <- y - (I-rho*W) %*% z/sigma^2 - exp(xb + z)/(1+exp(xb+z))
        dz
    }
    fit <- optim(rnorm(nn), logpost.Lap, gr=dlogpost.Lap, hessian = TRUE,
                 method = "BFGS", control = list(fnscale = -1) )
    QL <- -fit$hessian
    L <- chol(QL)

    t0 <- y %*% Z0 - sum(n.trial * log(1+exp(xb + Z0))) - Z0 %*% Q %*% Z0 / 2

    acc <- 0
    accept <- rep(0,N.Zy)

    if(scale.fixed){ ## fixed scale
        for(i in 1:N.Zy) {
            Zb <- Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(n.trial * log(1+exp(xb + Zb))) - Zb %*% Q %*% Zb / 2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <- acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z, accept))
    }
    else{ ## adaptive scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <- Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(n.trial * log(1+exp(xb + Zb))) - Zb %*% Q %*% Zb / 2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <- acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept,Scale.vec))
    }
}

########################################################################################
## Use MALA to sample from Z|y
postZ.mala.binom <- function(data, Z.start, psi, mcmc.pars, plots=FALSE){
#    for(i in 1:length(data)){
 #    tempobj=data[[i]]
  #   eval(parse(text=paste(names(data)[[i]],"= tempobj")))
 #}
#     for(i in 1:length(mcmc.pars)){
 #    tempobj=mcmc.pars[[i]]
  #   eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
# }
   	covX <- data$covX
	W <- data$W
 	y <- data$y
	n.trial <- data$n.trial
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

	Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    t0 <- y %*% Z0 - sum(n.trial * log(1 + exp(xb + Z0))) - Z0 %*% Q %*% Z0/2
    dlogpostZ0 <- as.numeric(y - Q %*% Z0 - n.trial * exp(xb + Z0)/(1 + exp(xb + Z0)))

    acc <- 0
    accept <- rep(0,N.Zy)

    if(scale.fixed){ ## Fixed Scale
        for(i in 1:N.Zy) {
            Scale2 <- Scale^2
            Zb <- Z0 + 0.5*Scale2*dlogpostZ0 + rnorm(nn, sd = Scale)
            dlogpostZb <- as.numeric(y - Q %*% Zb -
                                     n.trial * exp(xb + Zb)/(1 + exp(xb + Zb)))
            q10 <- -sum((Zb - Z0 - 0.5*Scale2*dlogpostZ0)^2) / (2*Scale2)
            q01 <- -sum((Z0 - Zb - 0.5*Scale2*dlogpostZb)^2) / (2*Scale2)
            t1 <- y %*% Zb - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb/2

            a <- as.numeric(exp(t1 + q01 - t0 - q10))
            r <- runif(1)
            if(r < a){
                Z0 <- Zb
                t0 <- t1
                dlogpostZ0 <- dlogpostZb
                acc <- acc + 1
            }
            accept[i] <- acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z, accept))
    }
    else{ ## Adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Scale2 <- Scale^2
            Zb <- Z0 + 0.5*Scale2*dlogpostZ0 + rnorm(nn, sd = Scale)
            dlogpostZb <- as.numeric(y - Q %*% Zb -
                                     n.trial * exp(xb + Zb)/(1 + exp(xb + Zb)))
            q10 <- -sum((Zb - Z0 - 0.5*Scale2*dlogpostZ0)^2) / (2*Scale2)
            q01 <- -sum((Z0 - Zb - 0.5*Scale2*dlogpostZb)^2) / (2*Scale2)
            t1 <- y %*% Zb - sum(n.trial * log(1 + exp(xb + Zb))) - Zb %*% Q %*% Zb/2

            a <- as.numeric(exp(t1 + q01 - t0 - q10))
            r <- runif(1)
            if(r < a){
                Z0 <- Zb
                t0 <- t1
                dlogpostZ0 <- dlogpostZb
                acc <- acc + 1
            }
            accept[i] <- acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.57)) # Adapative Scale so that the acceptence rate is 0.57
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z = Z, AC.r = accept, Scales = Scale.vec))
    }
}


#### Simulation from posterior from glm binomial with latent CAR
#### using different proposals
########################################################################################
postZ.rwmh.poisson <-  function(data, Z.start, psi, mcmc.pars, plots = FALSE){
   # for(i in 1:length(data)){
   #     tempobj=data[[i]]
   #     eval(parse(text=paste(names(data)[[i]],"= tempobj")))
   # }
   # for(i in 1:length(mcmc.pars)){
   #     tempobj=mcmc.pars[[i]]
   #     eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
   # }
 	W <- data$W
	covX <- data$covX
	y <- data$y
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

	Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    t0 <- sum(y * Z0) - sum(exp(xb + Z0)) - Z0 %*% Q %*% Z0 /2

    acc <- 0
    accept <- rep(0,N.Zy)

    if (scale.fixed){ ## fixed Scale
        for(i in 1:N.Zy) {
            Zb <- Z0 + rnorm(nn, sd = Scale)
            t1 <- sum(Zb * y) - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <-acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept))
    }
    else { ## adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <- Z0 + rnorm(nn, sd = Scale)
            t1 <- sum(Zb * y) - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <-acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept, Scale.vec))
    }
}


######################################################################################################
## Use RWHM with CAR proposal to sample from Z|y
postZ.CAR.poisson <-  function(data, Z.start, psi, mcmc.pars, plots = FALSE){
   # for(i in 1:length(data)){
   #     tempobj=data[[i]]
   #     eval(parse(text=paste(names(data)[[i]],"= tempobj")))
   # }
   # for(i in 1:length(mcmc.pars)){
   #     tempobj=mcmc.pars[[i]]
   #     eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
   # }
  	W <- data$W
	covX <- data$covX
	y <- data$y
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

	Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    L <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * nn)) ## upper triagle
    t0 <- y %*% Z0 - sum(exp(xb + Z0)) - Z0 %*% Q %*% Z0 / 2

    acc <- 0
    accept <- rep(0,N.Zy)

    ##if (scale.fix){
      if (scale.fixed){
      for(i in 1:N.Zy) { ## Fixed scale
            Zb <-Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc + 1
            }
            accept[i] <-acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept))
    }
    else{ ## Adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <-Z0 + solve(L, rnorm(nn, sd = scale))
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb /2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc + 1
            }
            accept[i] <-acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept, Scale.vec))
    }
}

#######################################################################################
## Use RWHM with Laplace approximation as proposal to sample from Z|y
postZ.Laplace.poisson  <-  function(data,  Z.start, psi, mcmc.pars, plots = FALSE){
    #for(i in 1:length(data)){
    #    tempobj=data[[i]]
    #    eval(parse(text=paste(names(data)[[i]],"= tempobj")))
    #}
    #for(i in 1:length(mcmc.pars)){
    #    tempobj=mcmc.pars[[i]]
    #    eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
    #}
  	W <- data$W
	covX <- data$covX
	y <- data$y
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

	Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma

    ## Laplace approximation to the posterior target
    logpost.Lap <- function(z){
        I <- as.matrix(diag.spam(1, nrow(W)))
        t1 <- sum((z+xb) %*% y)
        t2 <- -sum(exp(z+xb))
        t3 <- -0.5 * z %*% (I - rho*W) %*% z
        sum(t1,t2,t3)
    }
    dlogpost.Lap <- function(z){
        I <- as.matrix(diag.spam(1, nrow(W)))
        dz <- y - (I-rho*W) %*% z/sigma^2 - exp(xb + z)
        dz
    }
    fit <- optim(rnorm(nn), logpost.Lap, gr=dlogpost.Lap, hessian = TRUE,
                 method = "BFGS", control = list(fnscale = -1) )
    QL <- -fit$hessian
    L <- chol(QL)

    t0 <- y %*% Z0 - sum(exp(xb + Z0)) - Z0 %*% Q %*% Z0 / 2

    acc <- 0
    accept <- rep(0,N.Zy)

    if(scale.fixed){ ## fixed scale
        for(i in 1:N.Zy) {
            Zb <- Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb / 2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <- acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z, accept))
    }
    else{ ## adaptive scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Zb <- Z0 + solve(L, rnorm(nn, sd = Scale))
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb / 2

            a <- as.numeric(exp(t1 - t0))
            r <- runif(1)
            if (r < a){
                Z0 <- Zb
                t0 <- t1
                acc <- acc+1
            }
            accept[i] <- acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.24)) # Adapative Scale so that the acceptence rate is 0.24
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }

        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z,accept,Scale.vec))
    }
}

########################################################################################
## Use MALA to sample from Z|y
postZ.mala.poisson <- function(data, Z.start, psi, mcmc.pars, plots=FALSE){
    #for(i in 1:length(data)){
    #    tempobj=data[[i]]
    #    eval(parse(text=paste(names(data)[[i]],"= tempobj")))
    #}
    #for(i in 1:length(mcmc.pars)){
    #    tempobj=mcmc.pars[[i]]
    #    eval(parse(text=paste(names(mcmc.pars)[[i]],"= tempobj")))
    #}
  	W <- data$W
	covX <- data$covX
	y <- data$y
	c.ad <- mcmc.pars$c.ad
	thin <- mcmc.pars$thin
	burns <- mcmc.pars$burns
	N.Zy <- mcmc.pars$N.Zy
	Scale <- mcmc.pars$Scale
	scale.fixed <- mcmc.pars$scale.fixed
    Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

	Z0 <- Z.start
    rho <- psi[1]
    sigma <- psi[2]
    beta <- psi[-c(1,2)]
    c1.Scale <- c.ad[1]
    c2.Scale <- c.ad[2]

    nn <- nrow(W)
    Z <- matrix(0,ncol = nn, nrow =(N.Zy-burns)/thin)
    xb <- covX %*% beta
    I <- diag.spam(1, nn)
    Q <- (I - rho * W)/sigma
    t0 <- y %*% Z0 - sum(exp(xb + Z0)) - Z0 %*% Q %*% Z0/2
    dlogpostZ0 <- as.numeric(y - Q %*% Z0 - exp(xb + Z0))

    acc <- 0
    accept <- rep(0,N.Zy)

    if(scale.fixed){ ## Fixed Scale
        for(i in 1:N.Zy) {
            Scale2 <- Scale^2
            Zb <- Z0 + 0.5*Scale2*dlogpostZ0 + rnorm(nn, sd = Scale)
            dlogpostZb <- as.numeric(y - Q %*% Zb - exp(xb + Zb))
            q10 <- -sum((Zb - Z0 - 0.5*Scale2*dlogpostZ0)^2) / (2*Scale2)
            q01 <- -sum((Z0 - Zb - 0.5*Scale2*dlogpostZb)^2) / (2*Scale2)
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb/2

            a <- as.numeric(exp(t1 + q01 - t0 - q10))
            r <- runif(1)
            if(r < a){
                Z0 <- Zb
                t0 <- t1
                dlogpostZ0 <- dlogpostZb
                acc <- acc + 1
            }
            accept[i] <- acc/i
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z, accept))
    }
    else{ ## Adaptive Scale
        Scale.vec <- rep(0, N.Zy)
        for(i in 1:N.Zy) {
            Scale2 <- Scale^2
            Zb <- Z0 + 0.5*Scale2*dlogpostZ0 + rnorm(nn, sd = Scale)
            dlogpostZb <- as.numeric(y - Q %*% Zb - exp(xb + Zb))
            q10 <- -sum((Zb - Z0 - 0.5*Scale2*dlogpostZ0)^2) / (2*Scale2)
            q01 <- -sum((Z0 - Zb - 0.5*Scale2*dlogpostZb)^2) / (2*Scale2)
            t1 <- y %*% Zb - sum(exp(xb + Zb)) - Zb %*% Q %*% Zb/2

            a <- as.numeric(exp(t1 + q01 - t0 - q10))
            r <- runif(1)
            if(r < a){
                Z0 <- Zb
                t0 <- t1
                dlogpostZ0 <- dlogpostZb
                acc <- acc + 1
            }
            accept[i] <- acc/i
            Scale <- max(1/(nn^2*i), Scale + c1.Scale * i^(-c2.Scale) * (accept[i] - 0.57)) # Adapative Scale so that the acceptence rate is 0.57
            Scale.vec[i] <- Scale
            if(i > burns & (i-burns) %% thin == 0){
                Z[(i-burns)/thin, ] <- Z0
            }
        }
        ## The diagnostic plots
        if(plots){
            par(mfrow = c(2,2))
            plot(accept, type = "l")
            plot(Scale.vec, type = "l")

            acf.plot<-acf(Z[,1], plot = FALSE)
            plot(acf.plot$lag, acf.plot$acf, type = "l", xlab = "lag",
                 ylab = "autocorrelation", ylim = c(-0.2, 1))
            for(i in 2:nn){
                acf.plot <- acf(Z[, i], plot = FALSE)
                lines(acf.plot$lag, acf.plot$acf)
            }

            plot(Z[, 1], type = "l", xlab = "n.iter", ylab = "Z.i",
                 ylim = c(-10,10), col = 8)
            for(i in 2:5){
                lines(Z[,sample(1:nn,1)], col = i)
            }
        }
        return(list(Z = Z, AC.r = accept, Scales = Scale.vec))
    }
}
