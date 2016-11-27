#### Prepare Monte Carlo samples
mcl.prep.dCAR <- function(psi, n.samples, data){
    simY <- replicate(n.samples, CAR.simLM(psi, data))
    y <- data$data.vec$y
    t10 <- logunnorm(psi, data, y)
    t20 <- apply(simY, 2, logunnorm, pars=psi, data=data)
    return(list(simY = simY, t10=t10, t20=t20))
}

#### Monte Carlo log likelihood ratio
mcl.dCAR <- function(pars, data, simdata, rho.cons=c(-0.249, 0.249), Evar = FALSE){ # Evar has to turn off when doing optimization!
    data.vec <- data$data.vec
    y <- data.vec$y
    simY <- simdata$simY
    n.sim <- ncol(simY)
    t10 <- simdata$t10
    t20 <- simdata$t20

    t11 <- logunnorm(pars=pars,data=data, Y=y)
    t21 <- apply(simY, 2, logunnorm, pars=pars, data=data)

    s1 <- t11 - t10
    s2 <- t21 - t20

    # scale s2 in case of under/overflow
    s2.mean <- mean(s2)
    s2s <- s2 - s2.mean
    s2s.expm <- mean(exp(s2s))
    ##s2m <- s2s.expm * exp(s2.mean) # The normalising constant
    ls2m <- s2.mean + log(s2s.expm) # The log normalising constant

    if(is.null(rho.cons)){
    mc.lr <- as.numeric(s1-ls2m) # The log likelihood ratio
    v.lr <- var(exp(s2s)/s2s.expm)/n.sim # The estimated variance of the mc.lr # (eq 2.2.4)
    }else{
      if(pars[1] > rho.cons[1] & pars[1] < rho.cons[2]){
        mc.lr <- as.numeric(s1-ls2m) # The log likelihood ratio
        v.lr <- var(exp(s2s)/s2s.expm)/n.sim # The estimated variance of the mc.lr # (eq 2.2.4)
      }else{
        mc.lr <- -1e8
        v.lr <- 1e8
      }
    }

    if(Evar){
        return(c(mc.lr, v.lr))
    }else{
        return(mc.lr)
    }
}

#### Monte Carlo profile log likelihood ratio for rho
mcl.profile.dCAR <- function(rho, data, simdata, rho.cons = c(-0.249, 0.249), Evar = FALSE){
    data.vec <- data$data.vec
    y <- data.vec$y
    simY <- simdata$simY
    n.sim <- ncol(simY)
    t10 <- simdata$t10
    t20 <- simdata$t20

    pars <- c(rho, sigmabeta(rho, data))
    t11 <- logunnorm(pars=pars, data=data, Y=y)
    t21 <- apply(simY, 2, logunnorm, pars=pars, data=data)

    s1 <- t11 - t10
    s2 <- t21 - t20

    # scale s2 in case of under/overflow
    s2.mean <- mean(s2)
    s2s <- s2 - s2.mean
    s2s.expm <- mean(exp(s2s))
    ##s2m <- s2s.expm * exp(s2.mean) # The normalising constant
    ls2m <- s2.mean + log(s2s.expm) # The log normalising constant


    if(is.null(rho.cons)){
    mc.lr <- as.numeric(s1-ls2m) # The log likelihood ratio
    v.lr <- var(exp(s2s)/s2s.expm)/n.sim # The estimated variance of the mc.lr
    }else{
      if(rho> rho.cons[1] & rho < rho.cons[2]){
        mc.lr <- as.numeric(s1-ls2m) # The log likelihood ratio
        v.lr <- var(exp(s2s)/s2s.expm)/n.sim # The estimated variance of the mc.lr
      }
      else{
        mc.lr <- -1e8
        v.lr <- 1e8
      }
    }

    
    if(Evar){
         return(c(mc.lr, v.lr))
    }else{
        return(mc.lr)
    }
}

## Variance of the Monte Carlo MLE
vmle.dCAR <- function(MLE, data, simdata){
    par <- MLE$par
    ##n.sim <- ncol(simdata$simY)
    A <- mc.gradlik.dCAR(par, data, simdata, Evar = 2)$grad.var
    B <- MLE$hessian
    Binv <- solve(B)
    Binv %*% A %*% Binv # The matrix are usually small so just compute directly
}

#### Exact asympototic variance of t the MC maximum value
Avar.lik.dCAR <- function(pars, psi, data, n.samples, Log = TRUE){
    data.vec <-data$data.vec
    y <- data.vec$y
    n <- length(y)

    rho <- pars[1]
    sigma <- pars[2]

    rho.psi <- psi[1]
    sigma.psi <- psi[2]

    sigma.new <- sigma*sigma.psi/(2*sigma.psi - sigma)
    rho.new <- (2*rho*sigma.psi - rho.psi*sigma)/(2*sigma.psi - sigma)

    W <- data$W
    I <- diag(1,n)
    Q <- (I - rho*W)/sigma
    Q.det <- det(Q)

    Q.psi <-(I - rho.psi * W)/sigma.psi
    Q.psi.det <- det(Q.psi)

    EI.square <- Q.psi.det/Q.det

    if(abs(rho.new) >= 0.25){
        Q.new.det <- 0
    }else if(sigma.new > 0){
        Q.new <- (I - rho.new * W)/sigma.new
        Q.new.det <- det(Q.new)
    }else{
        Q.new.det <- 0}

    EI2 <- sqrt(Q.psi.det/Q.new.det)

    var.I <- (EI2 - EI.square)/n.samples
    if(Log){
        var.logI <- var.I/EI.square
        var.logI
    } else
    {var.I}

}

###############################################################################################
#### Other functions                                                                       ####
###############################################################################################
#### The variance of the MC maximum value
MCMLE.var <- function(pars, data, simdata){
    y <- data$data.vec$y
    simY <- simdata$simY
    d1 <- D.log.unorm(pars, data, y)
    d2 <- apply(simY, 2, D.log.unorm, pars = pars, data = data)
    t20 <- simdata$t20
    t21 <- apply(simY, 2, logunnorm, pars=pars, data=data)
    fZ <- (d1 - d2)*exp(t21-t20)
    var(t(fZ))
}

#### evaluate unnormalised density
logunnorm <- function(pars, data, Y){
    data.vec <- data$data.vec
    W <- data$W
    rho <- as.numeric(pars[1])
    sigma <- as.numeric(pars[2])
    n <- length(Y)
    I <- diag(1,n)
    Q <- (I - rho * W)/sigma

    if(length(pars) > 2){
        X <-  model.matrix(y~.,data=data.vec)
        beta <- pars[-c(1,2)]
        res <- Y-X %*% beta
        - t(res) %*% Q %*% res/2
    }else{
        - t(Y) %*% Q %*% Y/2
    }
}

#### derivative of unnormalised density
D.log.unorm <- function(pars, data, y){ # functional1 of the sample
    data.vec <- data$data.vec
    n <- length(y)
    W <- data$W
    I <- diag(1,n)
    rho <- as.numeric(pars[1])
    sigma <- as.numeric(pars[2])

    if(length(pars) > 2){
        beta <- pars[-c(1,2)]
        X <-  model.matrix(y~.,data=data.vec)
        Xb <- X %*% beta
        res <- y-Xb
        dsigma <- t(res) %*% ((I - rho*W) %*% res)/2/sigma^2
        drho <-  t(res) %*% (W %*% res)/2/sigma
        dbeta <- t(X) %*%(I - rho*W) %*% res/sigma
        return(c(as.numeric(drho), as.numeric(dsigma), as.numeric(dbeta)))
    }else{
        dsigma <- t(y) %*% ((I - rho*W) %*% y)/2/sigma^2
        drho <- t(y) %*% (W %*% y)/2/sigma
        return(c(as.numeric(drho), as.numeric(dsigma)))
    }
}

#### Gradient of the MC log likelihod
mc.gradlik.dCAR <- function(pars, data, simdata, Evar = 0){ # Evar has to turn off when doing optimization!
    data.vec <- data$data.vec
    y <- data.vec$y
    simY <- simdata$simY
    n.sim <- ncol(simY)
    ##t10 <- simdata$t10
    t20 <- simdata$t20

    W <- data$W
    rho <- as.numeric(pars[1])
    sigma <- as.numeric(pars[2])
    n <- length(y)
    I <- diag(1,n)
    Q <- (I - rho * W)/sigma

    log.unnorm <- function(Y){
    if(length(pars) > 2){
        X <-  model.matrix(y~.,data=data.vec)
        beta <- pars[-c(1,2)]
        res <- Y-X %*% beta
        - t(res) %*% Q %*% res/2
    }else{
        - t(Y) %*% Q %*% Y/2
    }
}

    D.log.unorm <- function(y){
        if(length(pars) > 2){
            beta <- pars[-c(1,2)]
            X <-  model.matrix(y~.,data=data.vec)
            Xb <- X %*% beta
            res <- y-Xb
            dsigma <- t(res) %*% ((I - rho*W) %*% res)/2/sigma^2
            drho <-  t(res) %*% (W %*% res)/2/sigma
            dbeta <- t(X) %*%(I - rho*W) %*% res/sigma
            return(c(as.numeric(drho), as.numeric(dsigma), as.numeric(dbeta)))
        }else{
            dsigma <- t(y) %*% ((I - rho*W) %*% y)/2/sigma^2
            drho <- t(y) %*% (W %*% y)/2/sigma
            return(c(as.numeric(drho), as.numeric(dsigma)))
        }
    }

    ##t11 <- log.unnorm(Y=y)
    t21 <- apply(simY, 2, log.unnorm)
    ##s1 <- t11 - t10
    s2 <- t21 - t20
    ## scale s2 in case of overflow
    s2.mean <- mean(s2)
    s2s <- s2 - s2.mean
    s2s.expm <- mean(exp(s2s))
    ##s2m <- s2s.expm * exp(s2.mean) # The normalising constant
    ##ls2m <- s2.mean + log(s2s.expm) # The log normalising constant

    g10 <- D.log.unorm(y)
    g20 <- apply(simY,2,D.log.unorm)
    g20w <- g20*exp(s2s)/s2s.expm

    grad.all <- g10 - rowMeans(g20w)

    if(Evar == 1){
        v.grad <- var(t(g20w))/n.sim
        return(list(grad = grad.all, grad.var = v.grad))
    }else if(Evar == 2){
      Ws <- exp(s2s)/s2s.expm
      grad.psi <- sapply(Ws, function(x) x*g10)
      v.grad <- var(t(g20w) - t(grad.psi))/n.sim
      return(list(grad = grad.all, grad.var = v.grad))
    }else{ return(grad.all)
      }
}


#### Hessian of the MC log likelihod
mc.Hesslik.dCAR <- function(pars, data, simdata){
    n.p <- length(pars)
    data.vec <- data$data.vec
    y.obs <- data.vec$y
    simY <- simdata$simY
    n.sim <- ncol(simY)
    ##t10 <- simdata$t10
    t20 <- simdata$t20

    W <- data$W
    rho <- as.numeric(pars[1])
    sigma <- as.numeric(pars[2])
    n <- length(y.obs)
    I <- diag(1,n)
    Q <- (I - rho * W)/sigma

    log.unnorm <- function(y){
        if(n.p > 2){
            X <-  model.matrix(y ~.,data=data.vec)
            beta <- pars[-c(1,2)]
            res <- y-X %*% beta
            - t(res) %*% Q %*% res/2
        }else{
            - t(y) %*% Q %*% y/2
        }
    }

    D.log.unorm <- function(y){
        if(n.p > 2){
            beta <- pars[-c(1,2)]
            X <-  model.matrix(y ~.,data=data.vec)
            Xb <- X %*% beta
            res <- y-Xb
            dsigma <- t(res) %*% ((I - rho*W) %*% res)/2/sigma^2
            drho <-  t(res) %*% (W %*% res)/2/sigma
            dbeta <- t(X) %*%(I - rho*W) %*% res/sigma
            return(c(as.numeric(drho), as.numeric(dsigma), as.numeric(dbeta)))
        }else{
            dsigma <- t(y) %*% ((I - rho*W) %*% y)/2/sigma^2
            drho <- t(y) %*% (W %*% y)/2/sigma
            return(c(as.numeric(drho), as.numeric(dsigma)))
        }
    }

    D2.log.unorm <- function(y){
        D2 <- matrix(0, n.p, n.p)
        if(n.p > 2){
            beta <- pars[-c(1,2)]
            X <-  model.matrix(y~.,data=data.vec)
            Xb <- X %*% beta
            res <- y-Xb

            d2rho <- 0
            d2sigma <- - t(res) %*% ((I - rho*W) %*% res)/sigma^3
            d2beta <- - t(X) %*%(I - rho*W) %*% X/sigma

            drhosigma <- -t(res) %*% (W %*% res)/2/sigma^2
            drhobeta <-  -t(X) %*% W %*% res/sigma
            dsigmabeta <-  -t(X) %*%(I - rho*W) %*% res/sigma^2


            D2[1,1] <- d2rho
            D2[2,2] <- d2sigma
            D2[3:n.p, 3:n.p] <- d2beta
            D2[1,2] <- D2[2,1] <- drhosigma
            D2[1,3:n.p] <- D2[3:n.p, 1] <- drhobeta
            D2[2,3:n.p] <- D2[3:n.p, 2] <- dsigmabeta

            return(D2)
        }else{
            d2rho <- 0
            d2sigma <- - t(y) %*% ((I - rho*W) %*% y)/sigma^3
            drhosigma <- -t(y) %*% (W %*% y)/2/sigma^2

            D2[1,1] <- d2rho
            D2[2,2] <- d2sigma
            D2[1,2] <- D2[2,1] <- drhosigma

            return(D2)
        }
    }

    ###### D2.lik = D2.log.unorm(y) - D2.log.C(theta)
    D2y <- D2.log.unorm(y.obs)
    #### D2.log.C(\theta) = A + B - D
    #### A = E(D2.log.unorm(Y)*W(Y))
    ## Calculate the weight: W(Y) = ci/C
    t21 <- apply(simY, 2, log.unnorm)
    s2 <- t21 - t20
    s2.mean <- mean(s2)
    s2s <- s2 - s2.mean
    s2s.expm <- mean(exp(s2s))
    WY <- exp(s2s)/s2s.expm
    ## Calculate A = mean(D2.log.unorm(Y)*WY)
    aa <- sapply(1:n.sim, function(x) D2.log.unorm(simY[,x])* WY[x])
    A <-  matrix(rowMeans(aa), n.p, n.p)
    ## Calculate B = mean(D.log.unorm(Y)^2 * WY)
    DY2 <- sapply(1:n.sim, function(x) WY[x] * crossprod(t(D.log.unorm(simY[,x]))))
    B <-  matrix(rowMeans(DY2 ), n.p, n.p)
    ## Calculate D = mean(D.log.unorm(Y) * WY)^2
    dd <- sapply(1:n.sim, function(x) D.log.unorm(simY[,x]) * WY[x])
    DYm <- rowMeans(dd)
    D <- crossprod(t(DYm))
    ## Finally
    D2.log.C <- A + B -D
    D2.lik <- D2y - D2.log.C
    return(D2.lik)
}
