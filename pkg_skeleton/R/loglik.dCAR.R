#####################################################################################
####                       Functions for direct CAR models                       ####
#####################################################################################
#### Exact likelihood evaluations for direct CAR models
#####################################################################################
#### exact log-likelihood function
loglik.dCAR <- function(pars,data){
    data.vec <-data$data.vec
    y <- data.vec$y
    n <- length(y)

    rho <- pars[1]
    sigma <- pars[2]

    if(length(pars) > 2){
    beta <- pars[-c(1,2)]
    X <- model.matrix(y~.,data=data.vec)
    Xb <- X %*% beta
    res <- y-Xb
}else{
    res <- y
}
    W <- data$W
    if(is.null(data$lambda)){
        lambda <-  eigen(W, symmetric = TRUE, only.values=TRUE)$values
    }else{
        lambda <- data$lambda
    }
    I <- diag(1,n)

    ll <- -n/2*(log(2*pi)+log(sigma)) + 1/2*sum(log(1 - rho*lambda)) -
          t(res) %*%((I - rho*W) %*% res)/(2 * sigma)
    as.numeric(ll)
}

#### The exact profile likelihood for rho
ploglik.dCAR <- function (rho, data){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    if(is.null(data$lambda)){
        lambda <-  eigen(W, symmetric = TRUE, only.values=TRUE)$values
    }else{
        lambda <- data$lambda
    }
    I <- diag(1,n)
    Q1 <- I - rho*W

    if(ncol(data.vec) > 1){
        X <- model.matrix(y~.,data=data.vec)
        beta <- solve(t(X) %*% Q1 %*% X) %*% t(X) %*% Q1 %*% y
        Xb <- X %*% beta
        res <- y - Xb
        sigma <- t(res) %*% Q1 %*% res/n
    }else{
        sigma <- t(y) %*% Q1 %*% y/n
    }
        ll <- -n/2*log(sigma) + 1/2*sum(log(1 - rho*lambda)) - n/2
        as.numeric(ll)
    }

#### Given rho get the mle of sigma and beta
sigmabeta <- function(rho, data){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    I <- diag(1,n)
    Q1 <- I - rho*W

    if(ncol(data.vec) > 1){
        X <- model.matrix(y~.,data=data.vec)
        beta <- solve(t(X) %*% Q1 %*% X) %*% t(X) %*% Q1 %*% y
        Xb <- X %*% beta
        res <- y - Xb
        sigma <- t(res) %*% Q1 %*% res/n
        return(c(sigma=as.numeric(sigma), beta=as.numeric(beta)))
    }else{
        sigma <- t(y) %*% Q1 %*% y/n
         return(as.numeric(sigma))
    }
}

#### Given rho get the mle of beta only
get.beta.lm <- function(rho, data){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    I <- diag(1,n)
    Q1 <- I - rho*W

    X <- model.matrix(y ~ .,data=data.vec)
    beta <- solve(t(X) %*% Q1 %*% X) %*% t(X) %*% Q1 %*% y
    beta

}

#### Maximum Pseudo-Likelihood estimator
mple.dCAR <- function(data, tol = 1e-6, rho0 = 0){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    I <- diag(1,n)
    diffrho = 1
    if(is.null(data$lambda)){
        lambda <-  eigen(W, symmetric = TRUE, only.values=TRUE)$values
    }else{
        lambda <- data$lambda
    }
    rr1 <- 1/min(lambda)
    rr2 <- 1/max(lambda)

    if(ncol(data.vec)>1){
        X <-  model.matrix(y~.,data=data.vec)
        while(diffrho >= tol){
            Q1 <- I-rho0*W
            beta0 <- solve(t(X) %*% Q1 %*% X) %*% (t(X) %*% Q1 %*% y)

            Xb <- X %*% beta0
            res <- y-Xb
            sigma0 <- as.numeric(t(res) %*% ((I - rho0*W) %*% res)/n)

            rho <- as.numeric(t(res) %*% (W %*% res)/(t(res) %*% (W %*% (W %*% res))))
            rho <- ifelse(abs(rho) >=0.25, sign(rho)*0.245, rho)
            diffrho <- abs(rho-rho0)
            rho0 <- rho
        }
            return(c(rho0, sigma0, beta0))
    }else{
        rho0 <- as.numeric(t(y) %*% (W %*% y)/(t(y) %*% (W %*% (W %*% y))))
        rho0 <- ifelse(rho0 < rr1, rr1*0.99, ifelse(rho0 > rr2, rr2*0.99, rho0))
        sigma0 <- as.numeric(t(y) %*% ((I - rho0*W) %*% y)/n)
        return(c(rho0, sigma0))
    }
}


#### gradient of the exact log-likelihood function
dloglik.dCAR <- function(pars,data){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    if(is.null(data$lambda)){
        lambda <-  eigen(W, symmetric = TRUE, only.values=TRUE)$values
    }else{
        lambda <- data$lambda
    }
    I <- diag(1,n)
    rho <- pars[1]
    sigma <- pars[2]

    if(length(pars) > 2){
        beta <- pars[-c(1,2)]
        X <-  model.matrix(y~.,data=data.vec)
        Xb <- X %*% beta
        res <- y-Xb
        dsigma <- -n/2/sigma + t(res) %*% ((I - rho*W) %*% res)/2/sigma^2
        drho <- -1/2*sum(lambda/(1 - rho*lambda)) + t(res) %*% (W %*% res)/2/sigma
        dbeta <- t(X) %*%(I - rho*W) %*% res/sigma
        return(c(as.numeric(drho), as.numeric(dsigma), as.numeric(dbeta)))
    }else{
        dsigma <- -n/2/sigma + t(y) %*% ((I - rho*W) %*% y)/2/sigma^2
        drho <- -1/2*sum(lambda/(1 - rho*lambda)) + t(y) %*% (W %*% y)/2/sigma
        return(c(as.numeric(drho), as.numeric(dsigma)))
    }
}

#### Hessian Matrix of the exact log-likelihood function
Hessian.dCAR <- function(pars,data){
    data.vec <- data$data.vec
    y <- data.vec$y
    n <- length(y)
    W <- data$W
    if(is.null(data$lambda)){
        lambda <-  eigen(W, symmetric = TRUE, only.values=TRUE)$values
    }else{
        lambda <- data$lambda
    }
    I <- diag(1,n)
    rho <- pars[1]
    sigma <- pars[2]

    if(length(pars) > 2){
        H.p <- length(pars)
        Hess <- matrix(0, nrow =H.p, ncol = H.p)
        beta <- pars[-c(1,2)]
        X <-  model.matrix(y~.,data=data.vec)
        Xb <- X %*% beta
        res <- y-Xb

        dsigma2 <- n/(2*sigma^2) - t(res) %*% ((I - rho*W) %*% res)/sigma^3
        drho2 <- -1/2*sum((lambda/(1 - rho*lambda))^2)
        dbeta2 <- -t(X)%*%(I - rho*W) %*% X/sigma

        drho.sig <- -t(res) %*% W %*% res/(2*sigma^2)
        drho.beta <- -t(X) %*% res/sigma
        dsig.beta <- -t(X) %*% (I - rho*W) %*% res/sigma^2

        Hess[1,1] <- drho2
        Hess[1,2] <- drho.sig
        Hess[2,1] <- drho.sig
        Hess[1, 3:H.p] <- drho.beta
        Hess[3:H.p, 1] <- drho.beta

        Hess[2,2] <- dsigma2
        Hess[2, 3:H.p] <- dsig.beta
        Hess[3:H.p, 2] <- dsig.beta

        Hess[3:H.p, 3:H.p] <- dbeta2
        return(Hess)
    }else{
        dsigma2 <- n/(2*sigma^2) - t(y) %*% ((I - rho*W) %*% y)/sigma^3
        drho2 <- -1/2*sum((lambda/(1 - rho*lambda))^2)
        drhosig <- -t(y) %*% W %*% y/(2*sigma^2)
        Hess <- matrix(c(drho2, drhosig, drhosig, dsigma2), nrow = 2, ncol = 2)
        return(Hess)
    }
}
