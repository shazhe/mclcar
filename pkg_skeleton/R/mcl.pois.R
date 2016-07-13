#####################################################################################
####                     Functions for glm Poisson models                        ####
#####################################################################################
## The joint density function of f(Y, Z) - poisson
logH.psi.pois <- function(pars, data, z){
    covX <- data$covX
    W <- data$W
    ##n.trial <- data$n.trial
    n.obs <- nrow(W)
    y <- data$y

    rho <- pars[1]
    sigma <- pars[2]
    beta <- as.vector(pars[-c(1,2)])

    xb <- covX %*% beta
    eta <- xb + z

    I <- diag(1, n.obs)
    Q <- (I - rho*W)/sigma
    Q1 <- as.spam(Q)
    L <- chol.spam(Q1, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n.obs))
    logdetQ.half <- sum(log(diag(L)))
    Const <- -n.obs/2*log(2*pi)+logdetQ.half

    tt <- Const - sum(sapply(y,lfactorial)) + crossprod(y, eta) -
        sum(exp(eta)) - 0.5*crossprod(z, crossprod(Q, z))
    return(tt)
}

########################################################################################
## Monte Carlo log-likelihood ratio
mcl.pois<- function(pars, mcdata, Evar = FALSE){
    rho <- pars[1]
    sigma <- pars[2]
    beta <- as.vector(pars[-c(1,2)])

    lHZy.psi <- mcdata$lHZy.psi
    Zy <- mcdata$Zy
    n.samples <- nrow(Zy)
    y <- mcdata$y
    W <- mcdata$W

    covX <- mcdata$covX
    ##n.trial <- mcdata$n.trial
    n.obs <- nrow(W)
    xb <- covX %*% beta ## possible speedup point with options
    I <- diag(1, n.obs)
    Q <- (I - rho*W)/sigma
    L <- chol(as.spam(Q), pivot = "MMD", memory = list(nnzcolindices = 6.25 * n.obs))
    logdetQ.half <- sum(log(diag(L)))
    Const <- -n.obs/2*log(2*pi)+logdetQ.half

    logH <- function(z){
        eta <- xb + z
        tt <- Const - sum(sapply(y,lfactorial)) + crossprod(y, eta) -
            sum(exp(eta)) - 0.5*crossprod(z, crossprod(Q, z))
        return(tt)
    }

    lHZy.pars <-apply(Zy, 1, logH)
    r.Zy <- lHZy.pars - lHZy.psi
    ## Scale with mean to avoid under/overflow in calculating the exponentials
    r.Zy.mean <- mean(r.Zy)
    r.Zy.scale <- r.Zy - r.Zy.mean
    r.exp.mean <- mean(exp(r.Zy.scale))

    mc.lr <-log(r.exp.mean) + r.Zy.mean


    if (Evar){
        ## Assuming the samples are independent after thinning
        ## So a very rough estimate using the sample variance
        var(exp(r.Zy))/mean(exp(r.Zy))^2/n.samples
        var(exp(r.Zy.scale))/(r.exp.mean^2)/n.samples
        v.lr <- var(exp(r.Zy.scale - log(r.exp.mean)))/n.samples

        return(c(mc.lr, v.lr))
    }
    else{
        return(mc.lr)
    }
}

########################################################################################
## Given rho and sigma, find beta
get.beta.pois <- function(beta0, rho.sig , mcdata){
    grad.beta <- function(beta){
        pars.b <- c(rho.sig, beta)
        grad.b <- mcl.grad.pois(pars.b, mcdata)[-c(1,2)]
        if (any(abs(beta0 - beta) > 2)){ # penalise when beta is too far away from the beta0
            return(sign(grad.b)*1e5)
        }else{
            return(grad.b)
        }
    }
    nleqslv(beta0, grad.beta,method = "Newton", control= list(ftol = 1, maxit = 2 ))
}


## The PROFILE Monte Carlo log-likelihood ratio
mcl.pois.profile <- function(pars, beta0, mcdata, Evar = FALSE){
    rho <- pars[1]
    sigma <- pars[2]
    covX <- mcdata$covX
   ## n.trial <- mcdata$n.trial

    lHZy.psi <- mcdata$lHZy.psi
    Zy <- mcdata$Zy
    y <- mcdata$y
    W <- mcdata$W
    n.obs <- nrow(W)
    n.samples <- nrow(Zy)

    grad.beta <- function(beta){
        pars.b <- c(rho,sigma, beta)
        grad.b <- mcl.grad.pois(pars.b,mcdata)[-c(1,2)]
        if (any(abs(beta0 - beta) > 2)){ # penalise when beta is too far away from the beta0
            return(sign(grad.b)*1e5)
        }else{
            return(grad.b)
        }
    }

    beta.solve <- nleqslv(beta0, grad.beta, control= list(ftol = 1, maxit = 3))
    beta <- beta.solve$x
    grad.beta <- beta.solve$fvec

    xb <- covX %*% beta
    I <- diag(1, n.obs)
    Q <- (I - rho*W)/sigma
    L <- chol(as.spam(Q), pivot = "MMD", memory = list(nnzcolindices = 6.25 * n.obs))
    logdetQ.half <- sum(log(diag(L)))
    Const <- -n.obs/2*log(2*pi)+logdetQ.half

    logH <- function(z){
        eta <- xb + z
        tt <- Const - sum(sapply(y,lfactorial)) + crossprod(y, eta) -
            sum(exp(eta)) - 0.5*crossprod(z, crossprod(Q, z))
        return(tt)
    }

    lHZy.pars <-apply(Zy, 1, logH)
    r.Zy <- lHZy.pars - lHZy.psi
    ## Scale with mean to avoid under/overflow in calculating the exponentials
    r.Zy.mean <- mean(r.Zy)
    r.Zy.scale <- r.Zy - r.Zy.mean
    r.exp.mean <- mean(exp(r.Zy.scale))

    mc.lr <- log(r.exp.mean) + r.Zy.mean
    if (Evar){
        ## Assuming the samples are independent after thinning
        ## So a very rough estimate using the sample variance
        v.lr <- var(exp(r.Zy.scale - log(r.exp.mean)))/n.samples
        return(c(mc.lr, v.lr, beta, grad.beta))
    }
    else{
        return(c(mc.lr, beta, grad.beta))
    }
}

#######################################################################################
## The gradient function for binomial model
mcl.grad.pois <- function(pars, mcdata, Evar = FALSE){
    rho <- pars[1]
    sigma <- pars[2]
    beta <- as.vector(pars[-c(1,2)])
    lHZy.psi <- mcdata$lHZy.psi
    Zy <- mcdata$Zy
    n.samples <- nrow(Zy)
    y <- mcdata$y
    W <- mcdata$W
    covX <- mcdata$covX
    ##n.trial <- mcdata$n.trial
    n.obs <- nrow(W)
    I <- diag(1, n.obs)
    Q <- (I - rho*W)/sigma

    xb <- covX %*% beta
    ew <- eigen(W, only.values = TRUE)$values
    logDet <- sum((1 - rho*ew))
    Const <- -n.obs/2*log(2*pi) + (logDet - n.obs * log(sigma))/2


    logH <- function(z){
        eta <- xb + z
        tt <- Const - sum(sapply(y,lfactorial)) + crossprod(y, eta) -
            sum(exp(eta)) - 0.5*crossprod(z, crossprod(Q, z))
        return(tt)
    }

    grad.logH <- function(z){
        eta <- z + xb

        g.rho <- 0.5 * z %*% W %*% z/sigma + 0.5*sum(-ew/(1-rho*ew))
        g.sigma <- -n.obs/(2*sigma) + (z %*% (I-rho*W) %*% z)/(2*sigma^2)
        g.beta <- colSums(as.numeric(y - exp(eta)) * covX)
        c(g.rho, g.sigma, g.beta)
    }

    lHZy.pars <-apply(Zy,1, logH)
    r.Zy <- lHZy.pars - lHZy.psi
    r.Zy.mean <- mean(r.Zy)        # Scale to deal with the overflow exp()
    r.Zy.scale <- r.Zy - r.Zy.mean
    W.Zy <- exp(r.Zy.scale)/mean(exp(r.Zy.scale))

    g.Zy.pars <- apply(Zy,1, grad.logH)
    grad.lli <- g.Zy.pars * W.Zy
    grad.ll <- rowMeans(grad.lli)

    if (Evar){
        ## Assuming the samples are independent after thinning
        ## So a very rough estimate using the sample variance
        v.lr <- var(t(grad.lli))/n.samples
        return(list(grad.lr = grad.ll, grad.var = v.lr))
    }
    else{
          return(grad.ll)
    }
}


########################################################################################
## The Hessian function
mcl.Hessian.pois <- function(pars, mcdata){
    rho <- pars[1]
    sigma <- pars[2]
    beta <- as.vector(pars[-c(1,2)])
    lHZy.psi <- mcdata$lHZy.psi
    Zy <- mcdata$Zy
    n.samples <- nrow(Zy)
    y <- mcdata$y
    W <- mcdata$W
    covX <- mcdata$covX
    ##n.trial <- mcdata$n.trial
    n.pars <- length(pars)
    n.obs <- nrow(W)
    I <- diag(1, n.obs)
    Q <- (I - rho*W)/sigma

    xb <- covX %*% beta
    ew <- eigen(W, only.values = TRUE)$values
    logDet <- sum((1 - rho*ew))
    Const <- -n.obs/2*log(2*pi) + (logDet - n.obs * log(sigma))/2

    logH <- function(z){
        eta <- xb + z
        tt <- Const - sum(sapply(y,lfactorial)) + crossprod(y, eta) -
            sum(exp(eta)) - 0.5*crossprod(z, crossprod(Q, z))
        return(tt)
    }

    grad.logH <- function(z){
        eta <- z + xb

        g.rho <- 0.5 * z %*% W %*% z/sigma + 0.5*sum(-ew/(1-rho*ew))
        g.sigma <- -0.5*n.obs/sigma + 0.5*(z %*% Q %*% z)/sigma
        g.beta <- colSums(as.numeric(y - exp(eta)) * covX)
        c(g.rho, g.sigma, g.beta)
    }

    Hessian.logH <- function(z){
        eta <- z + xb
        ##pz2 <- exp(eta)/(1 + exp(eta))^2
        H.mat <- matrix(0, n.pars, n.pars)

        H.rho2 <- -0.5*sum((ew/(1-rho*ew))^2)
        H.rho.sigma <- -0.5 * z %*% W %*% z/sigma^2
        H.sigma2 <-  -(z %*% Q %*% z)/(sigma^2) + 0.5*n.obs/sigma^2
        H.beta2 <- crossprod(as.numeric(exp(eta))*covX)
        H.mat[1,1] <- H.rho2
        H.mat[1,2] <- H.mat[2,1] <- H.rho.sigma
        H.mat[2,2] <- H.sigma2
        H.mat[-c(1,2), -c(1,2)] <- H.beta2
        H.mat
    }

    lHZy.pars <-apply(Zy,1, logH)
    r.Zy <- lHZy.pars - lHZy.psi
    r.Zy.mean <- mean(r.Zy)             # Scale to deal with the overflow exp()
    r.Zy.scale <- r.Zy - r.Zy.mean
    W.Zy <- exp(r.Zy.scale)/mean(exp(r.Zy.scale))

    ## -(grad(Z)*WZ)^2
    grad.wi <- sapply(1:n.samples, function(x) grad.logH(Zy[x,]) * W.Zy[x])
    grad.w <- rowMeans(grad.wi)
    H.part1 <- -crossprod(t(grad.w))
    ## mean(grad(Z)^2 * WZ)
    grad.2wi <- sapply(1:n.samples, function(x) crossprod(t(grad.logH(Zy[x,]))) * W.Zy[x])
    grad.2w <- rowMeans(grad.2wi)
    H.part2 <- matrix(grad.2w, n.pars, n.pars)
    ## mean(Hess(Z)*WZ)
    Hess.wi <- sapply(1:n.samples, function(x) Hessian.logH(Zy[x,]) * W.Zy[x])
    Hess.w <- rowMeans(Hess.wi)
    H.part3 <-  matrix(Hess.w, n.pars, n.pars)

    H <- H.part1 + H.part2 + H.part3
    H
}


############################################################
## Variance of the Monte Carlo MLE
vmle.pois <- function(MLE, mcdata){
    par <- MLE$estimate
    A <- mcl.grad.pois(par, mcdata, Evar = TRUE)$grad.var

    B <- MLE$hessian
    Binv <- solve(B)
    Binv %*% A %*% Binv # The matrix are usually small so just compute directly
}
